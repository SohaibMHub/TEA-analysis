# ===================================================================
#
#           FILE: tea_model.py
#           (The Backend Calculation Engine)
#
# This file contains the core process model class. It has no
# knowledge of the graphical interface and can be tested or
# used independently.
#
# ===================================================================

import numpy as np

class HydroPlantTEA:
    """
    A class to perform a Technoeconomic Analysis (TEA) for a 
    hydrometallurgical nickel recovery plant.
    It takes a dictionary of parameters and calculates all mass flows,
    energy consumption, CAPEX, OPEX, and key financial metrics.
    """
    def __init__(self, params):
        self.p = params
        self.r = {} # Results dictionary
        self._load_ore_composition()

    def _load_ore_composition(self):
        """
        Loads the fixed wt% composition of Serpentinised Peridotite.
        The data is normalized to ensure it sums to 100%. This is the
        single source of truth for the plant's feed.
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
        """Executes the calculation for each unit in the correct sequence."""
        self._calculate_unit_1_crushing()
        self._calculate_unit_2_mixing()
        self._calculate_unit_3_eleach()
        self._calculate_unit_4_centrifuge()
        self._calculate_unit_5_adsorption()
        self._calculate_unit_6_precipitation()
        self._calculate_final_economics()

    def _calculate_unit_1_crushing(self):
        """
        Calculates energy, power, size, and cost for the crushing/grinding circuit
        based on Bond's Law and engineering cost correlations.
        """
        # --- Constants and Conversions ---
        KG_PER_HR_TO_SHORT_TON_PER_HR = 1 / 907.185
        
        # --- Inputs from Parameters ---
        Wi = self.p['bond_work_index']
        P80 = self.p['p80_product_um']
        F80 = self.p['f80_feed_um']
        throughput_kghr = self.p['ore_flow_rate_kghr']
        
        # --- Energy & Power Calculation (Bond's Law) ---
        if P80 <= 0: P80 = 1 # Prevent division by zero
        
        # Specific Energy in kWh / short ton
        specific_energy_kwh_st = 10 * Wi * (1 / np.sqrt(P80) - 1 / np.sqrt(F80))
        throughput_stph = throughput_kghr * KG_PER_HR_TO_SHORT_TON_PER_HR
        theoretical_power_kw = specific_energy_kwh_st * throughput_stph
        electrical_power_kw = theoretical_power_kw / 0.90
        
        # --- CAPEX Calculation (Power Law Scaling) ---
        base_power = 1.6; base_capex = 130000; scaling_exponent = 0.6
        installed_capex_usd = base_capex * (electrical_power_kw / base_power) ** scaling_exponent if electrical_power_kw > 0 else 0

        # --- Conceptual Sizing (Heuristic) ---
        C_heuristic = 0.8 
        mill_diameter_m = (electrical_power_kw / (C_heuristic * 1.5))**(1/3.5) if electrical_power_kw > 0 else 0
        mill_length_m = mill_diameter_m * 1.5

        self.r['U1'] = {
            'capex': installed_capex_usd, 'power_kw': electrical_power_kw,
            'specific_energy_kwh_st': specific_energy_kwh_st,
            'mill_diameter_m': mill_diameter_m, 'mill_length_m': mill_length_m,
        }

    def _calculate_unit_2_mixing(self):
        """Calculates parameters for the slurry mixing tank."""
        self.r['U2'] = {'capex': self.p['capex_mixing'], 'power_kw': self.p['power_mixing_kw']}

    def _calculate_unit_3_eleach(self):
        """Calculates performance and cost for the electrochemical leaching cell."""
        self.r['ni_in_feed_kghr'] = self.p['ore_flow_rate_kghr'] * (self.ore_composition['Ni'] / 100)
        ca_in_ore_pct = self.ore_composition['CaO'] * (40.08 / 56.08)
        mg_in_ore_pct = self.ore_composition['MgO'] * (24.31 / 40.30)
        
        ca_in_kghr = self.p['ore_flow_rate_kghr'] * (ca_in_ore_pct / 100)
        mg_in_kghr = self.p['ore_flow_rate_kghr'] * (mg_in_ore_pct / 100)
        
        ni_leached_kghr = self.r['ni_in_feed_kghr'] * (self.p['leach_eff_ni_pct'] / 100)
        ca_leached_kghr = ca_in_kghr * (self.p['leach_eff_ca_pct'] / 100)
        mg_leached_kghr = mg_in_kghr * (self.p['leach_eff_mg_pct'] / 100)

        e_moles_ni = (ni_leached_kghr*1000 / 58.69) * 2
        e_moles_ca = (ca_leached_kghr*1000 / 40.08) * 2
        e_moles_mg = (mg_leached_kghr*1000 / 24.31) * 2
        total_e_moles_leaching = e_moles_ni + e_moles_ca + e_moles_mg
        total_e_moles = total_e_moles_leaching / (self.p['current_eff_leach_pct'] / 100)
        total_current_amps = total_e_moles * 96485 / 3600
        
        dc_power_kw = (total_current_amps * self.p['cell_voltage_v']) / 1000
        power_kw = dc_power_kw / 0.95
        
        electrode_area_m2 = total_current_amps / self.p['current_density_am2']
        capex_usd = (electrode_area_m2 * self.p['capex_cell_per_m2']) * self.p['inst_factor_cell']
        
        self.r['U3'] = {'capex': capex_usd, 'power_kw': power_kw}
        self.r['ni_leached_kghr'] = ni_leached_kghr

    def _calculate_unit_4_centrifuge(self):
        """Calculates solid-liquid separation performance and nickel loss."""
        self.r['U4'] = {'capex': self.p['capex_centrifuge_base'], 'power_kw': 1.4}
        solids_in_kghr = self.p['ore_flow_rate_kghr'] * (1 - (self.p['leach_eff_ni_pct'] + self.p['leach_eff_ca_pct'] + self.p['leach_eff_mg_pct'])/300)
        brine_kghr = self.p['ore_flow_rate_kghr'] * (100 / self.p['slurry_solids_pct'] - 1)
        liquid_in_cake_kghr = (solids_in_kghr / (1 - self.p['cake_moisture_pct']/100)) * (self.p['cake_moisture_pct']/100)
        fraction_liquid_lost = liquid_in_cake_kghr / brine_kghr if brine_kghr > 0 else 0
        
        self.r['ni_lost_to_cake_kghr'] = self.r['ni_leached_kghr'] * fraction_liquid_lost
        self.r['ni_to_adsorption_kghr'] = self.r['ni_leached_kghr'] - self.r['ni_lost_to_cake_kghr']

    def _calculate_unit_5_adsorption(self):
        """Calculates performance and cost for the adsorption unit."""
        adsorbent_vol_L = (self.r['ni_to_adsorption_kghr']*1000 / self.p['adsorbent_capacity_gL']) if self.p['adsorbent_capacity_gL'] > 0 else 0
        capex = self.p['capex_adsorption_base_per_L'] * adsorbent_vol_L + self.p['capex_adsorption_fixed']
        self.r['U5'] = {'capex': capex, 'power_kw': 0.25}
        self.r['ni_lost_to_barren_kghr'] = self.r['ni_to_adsorption_kghr'] * (1 - self.p['adsorption_eff_pct']/100)
        self.r['ni_to_precip_kghr'] = self.r['ni_to_adsorption_kghr'] - self.r['ni_lost_to_barren_kghr']

    def _calculate_unit_6_precipitation(self):
        """Calculates performance for the final product recovery unit."""
        self.r['U6'] = {'capex': self.p['capex_precip'], 'power_kw': self.p['power_precip_kw']}
        self.r['ni_lost_to_precip_kghr'] = self.r['ni_to_precip_kghr'] * (1 - self.p['precip_eff_pct']/100)
        self.r['ni_final_product_kghr'] = self.r['ni_to_precip_kghr'] - self.r['ni_lost_to_precip_kghr']

    def _calculate_final_economics(self):
        """Aggregates all costs and calculates final economic metrics."""
        units = ['U1', 'U2', 'U3', 'U4', 'U5', 'U6']
        total_capex = sum(self.r[u]['capex'] for u in units)
        total_power_kw = sum(self.r[u]['power_kw'] for u in units)
        self.r['total_capex'] = total_capex
        self.r['total_power_kw'] = total_power_kw

        cost_electricity_yr = total_power_kw * (self.p['electricity_cost_mwh']/1000) * self.p['op_hours_yr']
        cost_reagents_yr = self.p['reagent_cost_per_ton_ore'] * (self.p['ore_flow_rate_kghr']/1000) * self.p['op_hours_yr']
        cost_maintenance_yr = self.r['total_capex'] * self.p['maintenance_pct_capex']
        total_opex_yr = cost_electricity_yr + self.p['labor_cost_yr'] + cost_reagents_yr + cost_maintenance_yr
        self.r['total_opex_yr'] = total_opex_yr
        
        annual_ni_prod_kg = self.r['ni_final_product_kghr'] * self.p['op_hours_yr']
        self.r['annual_ni_prod_kg'] = annual_ni_prod_kg
        
        if annual_ni_prod_kg > 0:
            crf = (self.p['discount_rate'] * (1 + self.p['discount_rate'])**self.p['plant_life_yr']) / \
                  ((1 + self.p['discount_rate'])**self.p['plant_life_yr'] - 1)
            annualized_capex = self.r['total_capex'] * crf
            self.r['lcom_usd_kg'] = (annualized_capex + total_opex_yr) / annual_ni_prod_kg
        else:
            self.r['lcom_usd_kg'] = np.inf
            
        annual_revenue = annual_ni_prod_kg * self.p['ni_price_kg']
        self.r['annual_profit'] = annual_revenue - total_opex_yr

        self.r['capex_breakdown'] = {u: self.r[u]['capex'] for u in units}
        self.r['opex_breakdown'] = {'Electricity': cost_electricity_yr, 'Labor': self.p['labor_cost_yr'], 'Reagents': cost_reagents_yr, 'Maintenance': cost_maintenance_yr}
