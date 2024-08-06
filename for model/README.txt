Citation: Transgenerational virulence: Maternal pathogen exposure reduces offspring fitness
Authors: Kristina M. McIntire (1), Marcin K. Dziuba (1), Elizabeth Haywood (2), Miles L. Robertson (3), Megan Vaandrager (1), Emma Baird (1), Fiona Corcoran (1), Michael H. Cortez *(3), Meghan A. Duffy (1)
Contact: (1) Department of Ecology and Evolutionary Biology, University of Michigan, Ann Arbor, Michigan
(2) Department of Mathematics, Florida State University, Tallahassee, Florida
(3) Department of Biological Science, Florida State University, Tallahassee, Florida
        * Author responsible for writing code: mcortez@fsu.edu

Date: February 2023

-----------------------------------------------------------------------------------

ORGANIZATION

Folder SCIPModel contains files for SCIP model (Figures S2-S5)

Folder SCIDPModel contains files for SCIDP model (Figures 2, S5-S9)

-----------------------------------------------------------------------------------

CODE

Folder: SCIPModel
        - model_NXYP.m
                Matlab function for computing derivatives of frequency form of SCIP model.
        - FigureS2_TotalDensity_NXYP_script.m
                Simulates model and generates relationships between equilibrium total density and parameters of compromised individuals.
        - FigureS3_InfectionPrevalence_NXYP_script.m
                Simulates model and generates relationships between equilibrium infection prevalence and parameters of compromised individuals.
        - FigureS4_InfectedDensity_NXYP_script.m
                Simulates model and generates relationships between equilibrium infected density and parameters of compromised individuals.
        - FigureS5_BistabiltyExample_script.m
                Show example of bistability by simulating model with different initial conditions.

Folder: SCIDPModel
        - Calculations_NXYZP.mw
                Maple file that computes sensitivities for NXYZP model and computes equilibrium values for some sensitivities in Figures S6-S8; equilibrium values are saved as txt files.  
        - model_NXYZP.m
                Matlab function for computing derivatives of frequency form of SCIDP model.
        - FigureS6_TotalDensity_NXYP_script.m
                Simulates model and generates relationships between equilibrium total density and parameters of compromised and decimated individuals; some equilibrium values are imported from txt files generated in Maple.
        - FigureS7_InfectionPrevalence_NXYP_script.m
                Simulates model and generates relationships between equilibrium infection prevalence and parameters of compromised and decimated individuals; some equilibrium values are imported from txt files generated in Maple.
        - FigureS8_InfectedDensity_NXYP_script.m
                Simulates model and generates relationships between equilibrium infected density and parameters of compromised and decimated individuals; some equilibrium values are imported from txt files generated in Maple.
        - FigureS9_DentiferaPredictions_script.m
                Compare simulations of model parameterized to D. dentifera system with and without TGV.
        


 
