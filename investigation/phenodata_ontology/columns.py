"""Explicit, context-sensitive schemas for the 39 RPMM Phenodata CSVs.

Column terms describe attribute categories, not RDF predicates.  Values are
annotated separately.  A declared term ID is a candidate until the runner has
checked it against the pinned ontology catalog.  Original headers are never
rewritten, and identifiers/group codes stay local rather than becoming classes.
"""

from __future__ import annotations

import re
from copy import deepcopy


def normalize_header(header: str) -> str:
    """Normalize only for selecting a schema rule; never replace source text."""
    return re.sub(r"[^a-z0-9]+", "", header.lower())


def _spec(standard_name, semantic_role, value_kind, column_term_id="", unit_id="",
          reference_event="", status="accepted", reason="", rule_id=""):
    return {
        "standard_name": standard_name,
        "semantic_role": semantic_role,
        "column_term_id": column_term_id,
        "value_kind": value_kind,
        "unit_id": unit_id,
        "reference_event": reference_event,
        "status": status,
        "reason": reason,
        "rule_id": rule_id or "column:" + standard_name,
        "evidence": "",
    }


_RULES = {}


def _register(headers, *args, **kwargs):
    definition = _spec(*args, **kwargs)
    for header in headers:
        key = normalize_header(header)
        if key in _RULES and _RULES[key] != definition:
            raise ValueError("Conflicting header rule: " + header)
        _RULES[key] = definition


# Local identity and ordering are typed metadata, not ontology-class assertions.
_register(["gcid"], "gc_assay_id", "identifier", "identifier", status="local",
          reason="Local GC measurement identifier; do not infer biological replicate identity.")
_register(["lcid"], "lc_assay_id", "identifier", "identifier", status="local")
_register(["MSp"], "ms_assay_id", "identifier", "identifier", status="local")
_register(["FDDAp"], "fdda_assay_id", "identifier", "identifier", status="local")
_register(["Sample ID", "sampleID", "Numbering"], "sample_id", "identifier", "identifier", status="local")
_register(["sample name", "Sample name", "samplename", "sample No."],
          "sample_name", "identifier", "identifier", status="local",
          reason="Local name may encode several conditions; do not automatically split on punctuation.")
_register(["Order name", "Ordername", "ordername"], "source_order_name", "identifier", "identifier", status="local",
          reason="A source label, distinct from numeric acquisition runorder.")
_register(["MATLAB File name", "MATLABFilename", "MATLAB_Filename", "MATLABname"],
          "matlab_file_name", "file_reference", "identifier", status="local")
_register(["filename"], "raw_file_name", "file_reference", "identifier", status="local")
_register(["SIMCA-P obsID"], "simca_observation_id", "identifier", "identifier", status="local")
_register(["caroid"], "carotenoid_assay_id", "identifier", "identifier", status="local")
_register(["tube No.", "tube.No.", "tubeID", "tubeid", "tube name"],
          "tube_id", "identifier", "identifier", status="local",
          reason="Container label; some labels include source material descriptions.")
_register(["plant_id"], "source_plant_id", "identifier", "identifier", status="local")
_register(["No.", "num"], "source_row_number", "order", "integer", status="local")
_register(["runorder"], "acquisition_order", "order", "integer", status="local")
_register(["runorder_gc"], "gc_acquisition_order", "order", "integer", status="local")
_register(["runorder_lc"], "lc_acquisition_order", "order", "integer", status="local")
_register(["batch"], "analytical_batch", "group", "identifier", status="local")
_register(["proj"], "experiment_subgroup", "group", "identifier", status="local")
_register(["series"], "experiment_series", "group", "identifier", status="local")
_register(["class", "Class", "class in SIMCA", "class1", "class2", "mutantID"],
          "experimental_group_code", "group", "identifier", status="local",
          reason="Group labels are scoped to project/file/column; numeric class codes are not universal concepts.")

# Material properties and genetic identity.
_register(["species", "Species", "Speacies", "plant", "Plant", "sampletype"],
          "organism", "organism", "ontology_term", "COB:0000022",
          reason="A taxon value is resolved from the row; no project-wide taxon inheritance.")
_register(["speciesname"], "organism_common_name", "organism", "ontology_term", "COB:0000022",
          reason="Common name must be checked against the accompanying genus/epithet columns.")
_register(["binomialname"], "taxon_specific_epithet", "organism", "ontology_term", "COB:0000022",
          reason="Values are epithets, not full binomials; require the row's genus context.")
_register(["genotype", "Genotype", "geno", "geno1", "genotype1", "Genotype1",
           "genotype2", "Genotype2", "Genotype3"],
          "genotype_or_line", "genotype", "identifier", "EFO:0000513", status="local",
          reason="Preserve line/transgene/allele spelling; an external accession requires independent evidence.")
_register(["line", "Line", "Line name", "line number"], "plant_line", "genotype", "identifier", status="local",
          reason="Local plant-line accession; no assumed stock-center or gene database namespace.")
_register(["cultivar", "allium_cultivar"], "cultivar", "genotype", "identifier", status="local",
          reason="Cultivar names are kept distinct from genotype and taxon identity.")
_register(["genericback"], "genetic_background", "genotype", "identifier", status="local",
          reason="Source spelling genericback is preserved; MTB background label is a local identifier.")
_register(["gene loci"], "gene_locus", "identifier", "identifier", status="local",
          reason="Locus identifiers can link to a verified gene database; values such as control are not gene IDs.")
_register(["tissue", "Tissue", "fruitvsleaf"], "sampled_anatomical_structure", "organism_part", "ontology_term",
          "EFO:0000635", reason="Resolve plant structures with PO; reference/standard rows are handled separately.")
_register(["tissue1"], "anatomical_structure_level_1", "organism_part", "ontology_term", "EFO:0000635",
          reason="Top anatomical level; retain its link to tissue2 and tissue3 in the same row.")
_register(["tissue2"], "anatomical_structure_level_2", "organism_part", "ontology_term", "EFO:0000635",
          reason="Anatomical subtype inside tissue1; do not flatten the source hierarchy.")
_register(["tissue3"], "anatomical_structure_level_3", "organism_part", "ontology_term", "EFO:0000635",
          reason="Tissue/compartment inside tissue2; abbreviations need context.")
_register(["organ"], "anatomical_structure_or_sampling_position", "organism_part", "ontology_term", "EFO:0000635",
          status="needs_review", reason="Column mixes organ parts with composite leaf-position descriptions.")
_register(["sample_position_per_compound_leaf"], "position_relative_to_truss", "text", "text", status="local",
          reason="Ordinal and relative sampling position; retain original anatomical relationships.")
_register(["sampleposition_per_leaffoliole"], "leaflet_ordinal_position", "text", "text", status="local")
_register(["leaftype"], "relative_leaf_age", "phenotype", "ontology_term", "EFO:0000399",
          status="needs_review", reason="Old/middle are relative leaf-age categories, not defined absolute developmental stages.")
_register(["Color"], "organ_color", "phenotype", "ontology_term", "PATO:0000014")
_register(["Phenotype"], "organ_shape", "phenotype", "ontology_term", "PATO:0000052",
          reason="RPMM0020 values describe shape; affected organ must remain linked to the row context.")
_register(["visiblephenotype"], "visible_phenotype_present", "boolean", "boolean",
          reason="TRUE/FALSE indicate visible phenotype presence; not a phenotype class.")
_register(["pollination"], "pollination_performed", "boolean", "boolean",
          reason="TRUE/FALSE describes whether pollination was performed.")

# Quantities: generic mass is accepted without inventing what was weighed.
_register(["weight_mg", "weightmg", "weight (mg)"], "recorded_mass", "quantity", "quantity", "PATO:0000125",
          "UO:0000022", reason="Header explicitly specifies mg; sample/extract/injection-equivalent basis is not assumed.")
_register(["samplematerial_gFW"], "fresh_mass", "quantity", "quantity", "PATO:0000125", "UO:0000021",
          reason="Header specifies grams fresh weight.")
_register(["weight_mgFW_bulkedpistils2013"], "pooled_pistil_fresh_mass", "quantity", "quantity", "PATO:0000125",
          "UO:0000022", reason="Header specifies mg fresh weight of pooled pistils; source-year qualifier retained.")
_register(["plantweight_gc", "plantweight_lc"], "pooled_plant_mass", "quantity", "quantity", "PATO:0000125",
          status="needs_review", reason="Mass is indicated but no unit or fresh/dry basis is stated in this column.")
_register(["weightperplant_gc", "weightperplant_lc"], "mass_per_plant", "quantity", "quantity", "PATO:0000125",
          status="needs_review", reason="Per-plant mass is indicated but its unit/fresh-dry basis requires protocol evidence.")
_register(["injection volume (mg DW)"], "injected_dry_mass_equivalent", "quantity", "quantity", "PATO:0000125",
          "UO:0000022", status="needs_review",
          reason="Header calls the quantity volume but uses mg DW; mass-unit parsing is possible, physical meaning needs review.")
_register(["Area_cm2"], "sampled_surface_area", "quantity", "quantity", "PATO:0001323", "UO:0000081")
_register(["vialvolume_ml"], "vial_volume", "quantity", "quantity", "PATO:0000918", "UO:0000098")
_register(["DD_usedVolume_ul"], "standard_solution_volume", "quantity", "quantity", "PATO:0000918", "UO:0000101")
_register(["total volme for GC_ul"], "gc_extract_volume", "quantity", "quantity", "PATO:0000918", "UO:0000101")
_register(["conc_microM"], "standard_molar_concentration", "quantity", "quantity", "PATO:0000033", "UO:0000064")
_register(["conc"], "treatment_concentration", "quantity", "quantity", "PATO:0000033",
          reason="Unit must come from explicit cell suffix or same-row sample description; no default unit.")
_register(["DD_pmol"], "standard_substance_amount", "quantity", "quantity", unit_id="UO:0000042",
          reason="DD is a local abbreviation; amount unit pmol is explicit in header and values.")
_register(["alkaneconc_ng"], "alkane_standard_mass", "quantity", "quantity", "PATO:0000125", "UO:0000024",
          status="needs_review", reason="Header says concentration but provides ng without volume denominator; keep mass only.")
_register(["percentage_of_PEG"], "peg_percentage", "quantity", "quantity", unit_id="UO:0000187",
          status="needs_review", reason="Percent is explicit but mass/volume basis and PEG molecular-weight grade are not specified.")
_register(["mixtureratiofruit"], "fruit_mixture_fraction", "quantity", "quantity", unit_id="UO:0000187",
          status="needs_review", reason="Fruit and leaf fraction headers conflict with note; source must arbitrate before acceptance.")
_register(["mixtureratioleaf"], "leaf_mixture_fraction", "quantity", "quantity", unit_id="UO:0000187",
          status="needs_review", reason="Fruit and leaf fraction headers conflict with note; ratio basis also requires confirmation.")
_register(["nconc"], "relative_dilution_factor", "quantity", "quantity",
          reason="Unitless factor; values agree with concentration ratio strings in the adjacent column.")
_register(["concentration"], "dilution_ratio_description", "quantity", "quantity",
          reason="Values con: 1, con: 1/10, con: 1/100 describe dilution factors rather than molar concentrations.")
_register(["size"], "nominal_mass_class", "quantity", "quantity", "PATO:0000125", "UO:0000021",
          status="needs_review", reason="Values are nominal g classes/ranges including lessthan20g, not exact measured masses.")
_register(["incubationtemperature"], "incubation_temperature", "quantity", "quantity", "PATO:0000146",
          status="needs_review", reason="Temperature unit is absent; do not assume Celsius from the value 60.")
_register(["PPF_LED"], "photosynthetic_photon_flux", "quantity", "quantity",
          status="needs_review", reason="PPF lacks complete dimensional unit and area/time basis in the source header.")

_register(["num of samples combined per tube", "num_of_balkedsamples", "num_of_samples_bulked", "numofbulked",
           "bulkedsamples", "bulkgc"], "pooled_sample_count", "count", "integer",
          reason="Count of pooled source units, not a replicate identifier or biological replicate count.")
_register(["bulkedplants_gc"], "gc_pooled_plant_count", "count", "integer")
_register(["bulkedplants_lc"], "lc_pooled_plant_count", "count", "integer")
_register(["num of pistils"], "pooled_pistil_count", "count", "integer")

# Times retain their reference event; a unit alone cannot establish that event.
_register(["incubationtime_min"], "incubation_duration", "quantity", "quantity", "PATO:0001309", "UO:0000031",
          reference_event="incubation")
_register(["day_after_cutting"], "time_after_cutting", "quantity", "quantity", "PATO:0001309", "UO:0000034",
          reference_event="cutting", reason="Despite day in header, every source value explicitly uses weeks (w).")
_register(["floweringday"], "flowering_relative_time", "quantity", "quantity", "EFO:0000724", "UO:0000033",
          reference_event="flowering", reason="dbf/daf encode before/after flowering; sign must be retained.")
_register(["flowering_period"], "flowering_related_period", "quantity", "quantity", "PATO:0001309", "UO:0000034",
          status="needs_review", reason="Weeks are explicit but period direction and reference stage require protocol confirmation.")
_register(["time", "period", "timecourse", "day"], "recorded_elapsed_time", "quantity", "quantity", "PATO:0001309",
          status="needs_review", reason="Read units from each cell; the reference event is not established by this generic header.")
_register(["timeperiod"], "ordinal_sampling_stage", "text", "text", status="needs_review",
          reason="3rd/4th are ordinal labels; do not interpret them as a time duration.")
_register(["samplingdate"], "sampling_date", "date", "date",
          reason="YYYYMMDD source date can be represented as an ISO calendar date.")
_register(["samplingyear"], "sampling_year", "date", "date",
          reason="Year precision must be retained; do not invent month or day.")
_register(["preparedday"], "preparation_date", "date", "date", status="needs_review",
          reason="Six-digit YYMMDD value does not explicitly encode a century.")

# Experimental factors and protocol parameters.
_register(["treatment", "treatment1", "treatment2"], "treatment", "treatment", "ontology_term", "EFO:0000727",
          reason="Treatment concepts and their chemical agents must be represented separately.")
_register(["Ntreatment"], "nitrogen_treatment", "treatment", "ontology_term", "EFO:0000727")
_register(["Ptreatment"], "phosphorus_treatment", "treatment", "ontology_term", "EFO:0000727")
_register(["Inoculation", "infection2"], "inoculant_strain", "inoculant", "identifier", status="local",
          reason="Inoculant strain codes need a verified microbial taxon/strain accession; do not assign host genotype.")
_register(["infection1"], "inoculation_condition", "treatment", "ontology_term", "EFO:0000727")
_register(["innoculation"], "inoculation_performed", "boolean", "boolean",
          reason="Source typo retained; yes/no denotes performance of inoculation.")
_register(["growthcondition"], "growth_condition", "treatment", "ontology_term", "EFO:0000727")
_register(["nutrition"], "nutrition_regime_code", "treatment", "ontology_term", "EFO:0000727",
          status="needs_review", reason="NUT1000/NUT300 labels need medium/composition and dose definitions.")
_register(["light", "lightdark_period"], "illumination_condition", "treatment", "ontology_term", "EFO:0000727",
          reason="Light/dark/CL condition; do not infer intensity, photoperiod or duration without evidence.")
_register(["medium"], "incubation_medium", "treatment", "ontology_term",
          reason="Medium may be a mixture; preserve composition modifiers such as with sucrose.")
_register(["HS_device"], "headspace_sampling_method", "protocol_parameter", "text", status="local",
          reason="ITEX and SPME are sampling methods, not plant phenotypes.")
_register(["fiber"], "spme_fiber_specification", "protocol_parameter", "text", status="local",
          reason="Keep coating and thickness/product suffix; no unsupported ontology synonym collapse.")
_register(["cap"], "vial_cap_state", "protocol_parameter", "text", status="local")
_register(["exptype"], "experiment_type_code", "group", "identifier", status="local")
_register(["analysis"], "analytical_fraction", "protocol_parameter", "text", status="local",
          reason="cutin/wax describes the analysed fraction/method, not an organ name.")
_register(["new_oldanalysis"], "analysis_campaign", "group", "identifier", status="local")
_register(["type"], "sample_or_experiment_type_code", "group", "identifier", status="local")
_register(["combination"], "combined_factor_code", "group", "identifier", status="local",
          reason="Composite code summarizes adjacent tissue, nutrient and inoculant columns.")
_register(["what"], "anatomical_condition_code", "group", "identifier", status="local",
          reason="Abbreviations summarize adjacent organ/material descriptions; not universal concepts.")

# Descriptive fields are explicitly typed even where semantic decomposition needs review.
_register(["note", "note1", "note2", "notes"], "source_note", "text", "text", status="local")
_register(["client"], "sample_provider_label", "identifier", "identifier", status="local")
_register(["harvester"], "harvester_label", "identifier", "identifier", status="local")
_register(["annotation"], "gene_product_description", "text", "text", status="local",
          reason="Free text gene-product annotation; do not substitute a guessed GO term.")
_register(["Phenotype information "], "phenotype_description", "compound_description", "text", status="needs_review",
          reason="Free text combines line, generation, phenotype and control status; preserve entity-quality links.")
_register(["sample information", "sample information\n_Genotype_ Phenotype_ Tissue", "sampleinfo"],
          "combined_sample_description", "compound_description", "text", status="needs_review",
          reason="Several biological or experimental attributes are encoded together; use separate source columns when available.")
_register(["material"], "material_description", "compound_description", "text", status="needs_review",
          reason="Meaning varies by project: taxon, processed plant material or anatomical fraction.")
_register(["sample"], "sample_descriptor", "compound_description", "text", status="needs_review",
          reason="Meaning varies by project: sampling tool, organism or composite sample label.")
_register(["sepices"], "substrate_or_inoculant_taxon", "organism", "ontology_term", "COB:0000022",
          status="needs_review", reason="Column mixes soybean substrate and Bacillus inoculants; retain taxon role and strain separately.")
_register(["region"], "source_material_code", "text", "text", status="needs_review",
          reason="RPMM0062 value Soy is not a geographical region; exact source meaning needs confirmation.")


def classify_column(project: str, path: str, header: str) -> dict:
    """Return a fresh annotation specification for one original column.

    ``path`` may be relative to downloads/data or an absolute path.  Matching
    uses project plus path only for documented exceptions to the explicit
    header dictionary. Unknown headers always remain unresolved.
    """
    project = project.upper()
    key = normalize_header(header)
    result = deepcopy(_RULES.get(key))
    if result is None:
        result = _spec("unresolved_column", "unresolved", "unresolved", status="needs_review",
                       reason="No explicit schema rule for this header.", rule_id="column:unresolved")

    if project == "RPMM0037" and key == "genotype":
        result = _spec("fermentation_condition_code", "compound_description", "text", status="needs_review",
                       reason="Values such as 67_12 combine strain and time; this is not a plant-genotype field.",
                       rule_id="column:RPMM0037:genotype")
    elif project == "RPMM0003" and key == "sample":
        result = _spec("sampling_tool_or_reference_code", "protocol_parameter", "text", status="local",
                       reason="scalpel is a sampling tool; FD/STDs are reference-type codes, not organism names.",
                       rule_id="column:RPMM0003:sample")
    elif project == "RPMM0017" and key == "sample" and "valid01" in path.lower():
        result = _spec("organism_or_standard_type", "organism", "ontology_term", "COB:0000022",
                       reason="tomato/arabidopsis denote organisms; standard rows must be identified before taxon assignment.",
                       rule_id="column:RPMM0017:valid01:sample")
    elif project == "RPMM0004" and key == "material" and "/itex/" in ("/" + path.lower()):
        result = _spec("organism_genus", "organism", "ontology_term", "COB:0000022",
                       reason="allium provides genus context for the adjacent binomialname epithet.",
                       rule_id="column:RPMM0004:ITEX:material")
    elif project == "RPMM0036" and key == "material":
        result = _spec("sampled_material_layer", "organism_part", "ontology_term", "EFO:0000635",
                       reason="peel is the sampled anatomical fraction; retain parent Tissue=fruit.",
                       rule_id="column:RPMM0036:material")
    elif project == "RPMM0020" and key == "material":
        result = _spec("processed_source_material", "compound_description", "text", status="needs_review",
                       reason="Tomato Freeze dry combines taxon and material preparation; Plant supplies a separate taxon column.",
                       rule_id="column:RPMM0020:material")
    elif project == "RPMM0027" and key == "day" and "droughttreatment" in path.lower():
        result.update(standard_name="time_after_drought_treatment", unit_id="UO:0000033",
                      reference_event="", status="needs_review",
                      reason="DroughtTreatment path suggests treatment-relative days, but day zero requires protocol confirmation.",
                      rule_id="column:RPMM0027:DroughtTreatment:day")
    elif project == "RPMM0068" and key == "conc":
        result.update(standard_name="sodium_chloride_concentration",
                      reason="Resolve mM only when confirmed by the same row's sample name (0mM/200mM NaCl).",
                      rule_id="column:RPMM0068:conc")

    # This points to the exact source definition and does not claim external
    # term verification. The catalog check is a separate provenance layer.
    result["evidence"] = f"{path} :: header {header!r}; {result['reason']}"
    return result
