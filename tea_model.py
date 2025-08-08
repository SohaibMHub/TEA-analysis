# --- FILE 1: tea_model.py (The Backend Calculation Engine) ---
# This version is focused on Unit 1, using a fixed, complex ore composition.

import numpy as np

class HydroPlantTEA:
    """
    TEA model for a hydrometallurgical plant, starting with a fixed ore composition
    based on Serpentinised Peridotite.
    """
    def __init__(self, params):
        self.p = params
        self.r = {} # Results dictionary
        self._load_ore_composition()

    def _load_ore_composition(self):
        """
        Loads the fixed wt% composition of Serpentinised Peridotite.
        The data is normalized to sum to 100%.
        """
        # Data extracted from the user-provided image
        composition_data = {
            'SiO2': 41.65, 'TiO2': 0.034, 'Al2O3': 1.58, 'FeO': 7.83,
            'MnO': 0.125, 'MgO': 42.84, 'CaO': 0.58, 'Na2O': 0.0,
            'K2O': 0.0, 'P2O5': 0.016, 'LOI': 4.43
        }
        
        # Trace elements in ppm, converted to wt%
        trace_elements_ppm = {
            'Ni': 2547, 'Cr': 3154, 'V': 37, 'Sc': 10,
            'Cu': 18, 'Sr': 1, 'Y': 3, 'Zr': 6.2
        }
        for element, ppm in trace_elements_ppm.items():
            composition_data[element] = ppm / 10000.0 # Convert ppm to wt%

        # Normalize the composition to ensure it sums to 100%
        total_wt = sum(composition_data.values())
        self.ore_composition = {key: (value / total_wt) * 100 for key, value in composition_data.items()}
        self.r['ore_composition'] = self.ore_composition

    def run_full_analysis(self):
        """Executes the calculation for the defined units."""
        self._calculate_unit_1_crushing()
        # Placeholder for future units
        # self._calculate_unit_2_mixing()
        # ...etc...

    def _calculate_unit_1_crushing(self):
        """
        Calculates energy, power, and cost for the crushing/grinding circuit
        based on Bond's Law and engineering cost correlations.
        """
        # --- Energy & Power Calculation (Bond's Law) ---
        Wi = self.p['bond_work_index']
        P80 = self.p['p80_product_um']
        F80 = self.p['f80_feed_um']
        
        # Ensure P80 is not zero to avoid division errors
        if P80 <= 0: P80 = 1
        
        # Specific Energy Consumption in kWh / metric ton
        specific_energy_kwh_t = Wi * (10 / np.sqrt(P80) - 10 / np.sqrt(F80))
        
        # Throughput in metric tons / hour
        throughput_tph = self.p['ore_flow_rate_kghr'] / 1000
        
        # Theoretical power
        theoretical_power_kw = specific_energy_kwh_t * throughput_tph
        
        # Electrical power, assuming 90% motor efficiency
        electrical_power_kw = theoretical_power_kw / 0.90
        
        # --- CAPEX Calculation (Power Law Scaling) ---
        # Base case: 1.6 kW motor for a unit costing $130,000
        base_power = 1.6  # kW
        base_capex = 130000  # USD
        scaling_factor = 0.6 # Standard for this type of equipment
        
        # Use the scaling law, but prevent zero power from causing errors
        if electrical_power_kw > 0:
            installed_capex_usd = base_capex * (electrical_power_kw / base_power) ** scaling_factor
        else:
            installed_capex_usd = 0

        # Store all results for this unit
        self.r['U1'] = {
            'specific_energy_kwh_t': specific_energy_kwh_t,
            'electrical_power_kw': electrical_power_kw,
            'installed_capex_usd': installed_capex_usd
        }
