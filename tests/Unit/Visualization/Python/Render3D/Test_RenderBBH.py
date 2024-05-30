#!/usr/bin/env python

# Distributed under the MIT License.
# See LICENSE.txt for details.

import os
import shutil
import unittest

import numpy as np
from click.testing import CliRunner

from spectre.Informer import unit_test_build_path, unit_test_src_path
from spectre.Visualization.RenderBBH import render_bbh_command


class TestRenderBBH(unittest.TestCase):
    def setUp(self):
        self.vol_test_data = os.path.join(
            unit_test_src_path(), "Visualization/Python/VolTestData.xmf"
        )
        self.test_dir = os.path.join(
            unit_test_build_path(), "Visualization/Render3D/BBH"
        )
        self.output_file = os.path.join(self.test_dir, "test_clip.png")
        shutil.rmtree(self.test_dir, ignore_errors=True)
        os.makedirs(self.test_dir, exist_ok=True)

    def test_render_bbh(self):
        runner = CliRunner()
        result = runner.invoke(
            render_bbh_command,
            ["-o", self.test_dir, "-v", self.vol_test_data, "-c", "0"],
            catch_exceptions=False,
        )
        self.assertEqual(result.exit_code, 0, msg=result.output)
        self.assertTrue(os.path.exists(self.output_file), msg=result.output)


if __name__ == "__main__":
    unittest.main(verbosity=2)
