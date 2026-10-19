"""Physical invariants, independent reference checks and inverse recovery."""
from pathlib import Path
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"src"))
from moscap import Device, EPS_SI, reference_bias, reference_charge, solve_bias, sweep, system
from parameter_extraction import fit_cv


class MOSPhysicsTests(unittest.TestCase):
    def test_flat_band_and_bulk_neutrality(self):
        d = Device()
        s = solve_bias(d, d.vfb)
        np.testing.assert_allclose(s["psi"], 0, atol=1e-14)
        np.testing.assert_allclose(s["rho"], 0, atol=1e-14)
        self.assertAlmostEqual((d.p0-d.n0)/d.na, 1)
        csi = EPS_SI/d.ld*np.sqrt(1+d.eta)
        expected = d.cox*csi/(d.cox+csi)
        self.assertLess(abs(s["Cqs_F_m2"]/expected-1), 2e-4)

    def test_independent_first_integral_reference(self):
        d = Device()
        for v in (-.25, .4, 1.05):
            s, r = solve_bias(d, v), reference_bias(d, v)
            self.assertLess(abs(s["psi_s_V"]-r["psi_s_V"]), 2e-5)
            self.assertLess(abs(s["Cqs_F_m2"]-r["Cqs_F_m2"])/d.cox, 3e-4)

    def test_charge_conservation_and_neutral_far_boundary(self):
        d = Device()
        for v in (-.25, .4, 1.05):
            s = solve_bias(d, v)
            self.assertLess(abs(s["charge_balance_C_m2"]), 1e-10)
            self.assertLess(abs(s["bulk_displacement_C_m2"]), 1e-10)
            self.assertLess(s["residual_inf"], 2e-10)
            self.assertEqual(s["psi"][-1], 0)

    def test_surface_regimes(self):
        d = Device()
        a, dep, inv = [solve_bias(d, v) for v in (-.25, .4, 1.05)]
        self.assertGreater(a["p"][0], d.p0)
        self.assertLess(dep["p"][0], d.p0)
        self.assertLess(dep["psi_s_V"], 2*d.phi_f)
        self.assertGreater(inv["psi_s_V"], 2*d.phi_f)
        self.assertGreater(inv["n"][0], d.na)

    def test_mesh_refinement(self):
        d = Device()
        r = reference_bias(d, 1.05)
        errors = [abs(solve_bias(d, 1.05, n)["psi_s_V"]-r["psi_s_V"]) for n in (101, 201, 401)]
        self.assertTrue(3.7 < errors[0]/errors[1] < 4.3)
        self.assertTrue(3.7 < errors[1]/errors[2] < 4.3)

    def test_capacitance_is_charge_derivative(self):
        d = Device()
        step = 1e-4
        for v in (0., .4, 1.05):
            s = solve_bias(d, v)
            derivative = (solve_bias(d, v+step)["Qg_C_m2"]-solve_bias(d, v-step)["Qg_C_m2"])/(2*step)
            self.assertLess(abs(derivative-s["Cqs_F_m2"])/d.cox, 2e-6)

    def test_reference_derivative_at_and_near_flatband(self):
        d = Device()
        for psi in (-.1, -1e-8, 0., 1e-8, .5, .8):
            step = 1e-7
            numerical = -(reference_charge(d, psi+step)[0]-reference_charge(d, psi-step)[0])/(2*step)
            expected = reference_charge(d, psi)[1]
            self.assertLess(abs(numerical/expected-1), 1e-6)

    def test_frozen_electron_hf_limit(self):
        d = Device()
        for v in (-.25, .4, 1.05):
            s = solve_bias(d, v)
            self.assertGreater(s["Chf_frozen_F_m2"], 0)
            self.assertLessEqual(s["Chf_frozen_F_m2"], s["Cqs_F_m2"]+1e-14)
        inv = solve_bias(d, 1.05)
        self.assertLess(inv["Chf_frozen_F_m2"]/d.cox, .15)
        self.assertGreater(inv["Cqs_F_m2"]/d.cox, .6)

    def test_fixed_charge_shifts_flatband(self):
        base, charged = Device(), Device(qox_cm2=1e11)
        self.assertAlmostEqual(charged.vfb, -charged.qox/charged.cox)
        s0, s1 = solve_bias(base, .4), solve_bias(charged, .4+charged.vfb)
        np.testing.assert_allclose(s0["psi"], s1["psi"], atol=1e-10)
        self.assertAlmostEqual(s1["Qg_C_m2"]-s0["Qg_C_m2"], -charged.qox, places=12)
        self.assertLess(abs(s1["charge_balance_C_m2"]), 1e-10)

    def test_analytic_jacobian_directional_derivative(self):
        d = Device()
        s = solve_bias(d, .4, 101)
        u, y = s["u"], s["x"]/d.ld
        f, band = system(d, y, u, .4)
        vector = np.sin(np.linspace(0, 2, len(u)))
        product = band[1]*vector
        product[:-1] += band[0, 1:]*vector[1:]
        product[1:] += band[2, :-1]*vector[:-1]
        step = 1e-5
        numeric = (system(d, y, u+step*vector, .4)[0]-system(d, y, u-step*vector, .4)[0])/(2*step)
        np.testing.assert_allclose(numeric, product, rtol=3e-6, atol=3e-7)

    def test_continuation_and_direct_solve_agree(self):
        d = Device()
        result = sweep(d, np.array([-.25, 0., .4, 1.05]), 201)[-1]
        direct = solve_bias(d, 1.05, 201)
        np.testing.assert_allclose(result["psi"], direct["psi"], atol=1e-8)

    def test_noise_free_inverse_recovery(self):
        d = Device(na_cm3=1.2e16, tox_nm=12.)
        voltage = np.linspace(-.3, 1.1, 17)
        synthetic = np.array([reference_bias(d, v)["Cqs_F_m2"] for v in voltage])
        result = fit_cv(voltage, synthetic, .004*d.cox, d, nodes=201)
        self.assertLess(abs(result["device"].tox_nm/d.tox_nm-1), .01)
        self.assertLess(abs(result["device"].na_cm3/d.na_cm3-1), .01)

    def test_input_guards(self):
        for kwargs in ({"na_cm3": -1}, {"tox_nm": 0}, {"temperature": float("nan")}):
            with self.assertRaises(ValueError):
                Device(**kwargs)
        with self.assertRaises(ValueError):
            solve_bias(Device(), float("nan"))
        with self.assertRaises(ValueError):
            solve_bias(Device(), .4, nodes=10)


if __name__ == "__main__":
    unittest.main()
