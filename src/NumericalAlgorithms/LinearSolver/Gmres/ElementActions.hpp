// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <tuple>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/DenseVector.hpp"
#include "NumericalAlgorithms/LinearSolver/Gmres/ResidualMonitorActions.hpp"
#include "NumericalAlgorithms/LinearSolver/InnerProduct.hpp"
#include "NumericalAlgorithms/LinearSolver/Tags.hpp"
#include "Parallel/ConstGlobalCache.hpp"
#include "Parallel/Info.hpp"
#include "Parallel/Invoke.hpp"
#include "Parallel/Reduction.hpp"
#include "Utilities/Functional.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Requires.hpp"

/// \cond
namespace tuples {
template <typename...>
class TaggedTuple;
}  // namespace tuples
namespace LinearSolver {
namespace gmres_detail {
template <typename Metavariables>
struct ResidualMonitor;
}  // namespace gmres_detail
}  // namespace LinearSolver
/// \endcond

namespace LinearSolver {
namespace gmres_detail {

struct PrepareStep {
  template <typename DbTagsList, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ActionList,
            typename ParallelComponent>
  static std::tuple<db::DataBox<DbTagsList>&&> apply(
      db::DataBox<DbTagsList>& box,
      const tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      const Parallel::ConstGlobalCache<Metavariables>& /*cache*/,
      const ArrayIndex& /*array_index*/, const ActionList /*meta*/,
      const ParallelComponent* const /*meta*/) noexcept {
    using orthogonalization_iteration_id_tag =
        db::add_tag_prefix<LinearSolver::Tags::Orthogonalization,
                           LinearSolver::Tags::IterationId>;

    db::mutate<LinearSolver::Tags::IterationId,
               orthogonalization_iteration_id_tag>(
        make_not_null(&box),
        [](const gsl::not_null<db::item_type<LinearSolver::Tags::IterationId>*>
               iteration_id,
           const gsl::not_null<
               db::item_type<orthogonalization_iteration_id_tag>*>
               orthogonalization_iteration_id,
           const db::item_type<::Tags::Next<LinearSolver::Tags::IterationId>>&
               next_iteration_id) noexcept {
          *iteration_id = next_iteration_id;
          *orthogonalization_iteration_id = 0;
        },
        get<::Tags::Next<LinearSolver::Tags::IterationId>>(box));
    return {std::move(box)};
  }
};

struct PerformStep {
  template <typename DataBox, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ActionList,
            typename ParallelComponent>
  static std::tuple<DataBox&&, bool> apply(
      DataBox& box, const tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      Parallel::ConstGlobalCache<Metavariables>& cache,
      const ArrayIndex& array_index,
      // NOLINTNEXTLINE(readability-avoid-const-params-in-decls)
      const ActionList /*meta*/,
      // NOLINTNEXTLINE(readability-avoid-const-params-in-decls)
      const ParallelComponent* const /*meta*/) noexcept {
    using fields_tag = typename Metavariables::system::fields_tag;
    using operand_tag =
        db::add_tag_prefix<LinearSolver::Tags::Operand, fields_tag>;
    using operator_tag =
        db::add_tag_prefix<LinearSolver::Tags::OperatorAppliedTo, operand_tag>;
    using basis_history_tag =
        LinearSolver::Tags::KrylovSubspaceBasis<fields_tag>;

    ASSERT(get<LinearSolver::Tags::IterationId>(box) !=
               std::numeric_limits<size_t>::max(),
           "Linear solve iteration ID is at initial state. Did you forget to "
           "invoke 'PrepareStep'?");

    db::mutate<operand_tag>(
        make_not_null(&box),
        [](const gsl::not_null<db::item_type<operand_tag>*> operand,
           const db::item_type<operator_tag>& operator_action) noexcept {
          *operand = db::item_type<operand_tag>(operator_action);
        },
        get<operator_tag>(box));

    Parallel::contribute_to_reduction<
        StoreOrthogonalization<ParallelComponent>>(
        Parallel::ReductionData<
            Parallel::ReductionDatum<double, funcl::Plus<>>>{inner_product(
            get<basis_history_tag>(box)[0], get<operand_tag>(box))},
        Parallel::get_parallel_component<ParallelComponent>(cache)[array_index],
        Parallel::get_parallel_component<ResidualMonitor<Metavariables>>(
            cache));

    // Terminate algorithm for now. The `ResidualMonitor` will receive the
    // reduction that is performed above and then broadcast to the following
    // action, which is responsible for restarting the algorithm.
    return {std::move(box), true};
  }
};

struct OrthogonalizeOperand {
  template <
      typename ParallelComponent, typename DataBox, typename Metavariables,
      typename ArrayIndex,
      Requires<db::tag_is_retrievable_v<
                   typename Metavariables::system::fields_tag, DataBox> and
               db::tag_is_retrievable_v<LinearSolver::Tags::HasConverged,
                                        DataBox>> = nullptr>
  static void apply(DataBox& box,
                    Parallel::ConstGlobalCache<Metavariables>& cache,
                    const ArrayIndex& array_index,
                    const double orthogonalization) noexcept {
    using fields_tag = typename Metavariables::system::fields_tag;
    using operand_tag =
        db::add_tag_prefix<LinearSolver::Tags::Operand, fields_tag>;
    using orthogonalization_iteration_id_tag =
        db::add_tag_prefix<LinearSolver::Tags::Orthogonalization,
                           LinearSolver::Tags::IterationId>;
    using basis_history_tag =
        LinearSolver::Tags::KrylovSubspaceBasis<fields_tag>;

    db::mutate<operand_tag, orthogonalization_iteration_id_tag>(
        make_not_null(&box),
        [orthogonalization](
            const gsl::not_null<db::item_type<operand_tag>*> operand,
            const gsl::not_null<
                db::item_type<orthogonalization_iteration_id_tag>*>
                orthogonalization_iteration_id,
            const db::item_type<basis_history_tag>& basis_history) noexcept {
          *operand -= orthogonalization *
                      gsl::at(basis_history, *orthogonalization_iteration_id);
          (*orthogonalization_iteration_id)++;
        },
        get<basis_history_tag>(box));

    const auto& next_orthogonalization_iteration_id =
        get<orthogonalization_iteration_id_tag>(box);
    const auto& iteration_id = get<LinearSolver::Tags::IterationId>(box);

    if (next_orthogonalization_iteration_id <= iteration_id) {
      Parallel::contribute_to_reduction<
          StoreOrthogonalization<ParallelComponent>>(
          Parallel::ReductionData<
              Parallel::ReductionDatum<double, funcl::Plus<>>>{
              inner_product(gsl::at(get<basis_history_tag>(box),
                                    next_orthogonalization_iteration_id),
                            get<operand_tag>(box))},
          Parallel::get_parallel_component<ParallelComponent>(
              cache)[array_index],
          Parallel::get_parallel_component<ResidualMonitor<Metavariables>>(
              cache));
    } else {
      Parallel::contribute_to_reduction<
          StoreFinalOrthogonalization<ParallelComponent>>(
          Parallel::ReductionData<
              Parallel::ReductionDatum<double, funcl::Plus<>>>{
              inner_product(get<operand_tag>(box), get<operand_tag>(box))},
          Parallel::get_parallel_component<ParallelComponent>(
              cache)[array_index],
          Parallel::get_parallel_component<ResidualMonitor<Metavariables>>(
              cache));
    }
  }
};

struct NormalizeOperandAndUpdateField {
  template <
      typename ParallelComponent, typename DataBox, typename Metavariables,
      typename ArrayIndex,
      Requires<db::tag_is_retrievable_v<
                   typename Metavariables::system::fields_tag, DataBox> and
               db::tag_is_retrievable_v<LinearSolver::Tags::HasConverged,
                                        DataBox>> = nullptr>
  static void apply(DataBox& box,
                    Parallel::ConstGlobalCache<Metavariables>& cache,
                    const ArrayIndex& array_index, const double normalization,
                    const DenseVector<double>& minres,
                    const db::item_type<LinearSolver::Tags::HasConverged>&
                        has_converged) noexcept {
    using fields_tag = typename Metavariables::system::fields_tag;
    using initial_fields_tag =
        db::add_tag_prefix<LinearSolver::Tags::Initial, fields_tag>;
    using operand_tag =
        db::add_tag_prefix<LinearSolver::Tags::Operand, fields_tag>;
    using basis_history_tag =
        LinearSolver::Tags::KrylovSubspaceBasis<fields_tag>;

    db::mutate<operand_tag, basis_history_tag, fields_tag,
               LinearSolver::Tags::HasConverged>(
        make_not_null(&box),
        [
          normalization, &minres, &has_converged
        ](const gsl::not_null<db::item_type<operand_tag>*> operand,
          const gsl::not_null<db::item_type<basis_history_tag>*> basis_history,
          const gsl::not_null<db::item_type<fields_tag>*> field,
          const gsl::not_null<db::item_type<LinearSolver::Tags::HasConverged>*>
              local_has_converged,
          const db::item_type<initial_fields_tag>& initial_field) noexcept {
          *operand /= normalization;
          basis_history->push_back(*operand);
          *field = initial_field;
          for (size_t i = 0; i < minres.size(); i++) {
            *field += minres[i] * gsl::at(*basis_history, i);
          }
          *local_has_converged = has_converged;
        },
        get<initial_fields_tag>(box));

    // Proceed with algorithm
    Parallel::get_parallel_component<ParallelComponent>(cache)[array_index]
        .perform_algorithm(true);
  }
};

}  // namespace gmres_detail
}  // namespace LinearSolver
