"""Conservative, context-specific annotations for RPMM Phenodata cell values.

This module never changes the source value.  Identifier and quantity cells have
types rather than fabricated ontology class identifiers.  Term identifiers are
checked against the pinned ontology catalogue by the caller before publication.
"""

from datetime import datetime
from decimal import Decimal, InvalidOperation
from functools import lru_cache
import json
from pathlib import Path
import re


SPECIES = {
    "arabidopsis": ("NCBITaxon:3702", "Arabidopsis thaliana"),
    "arabidopsis thaliana": ("NCBITaxon:3702", "Arabidopsis thaliana"),
    "tomato": ("NCBITaxon:4081", "Solanum lycopersicum"),
    "solanum lycopersicum": ("NCBITaxon:4081", "Solanum lycopersicum"),
    "rice": ("NCBITaxon:4530", "Oryza sativa"),
    "oryza sativa": ("NCBITaxon:4530", "Oryza sativa"),
    "lettuce": ("NCBITaxon:4236", "Lactuca sativa"),
    "cassava": ("NCBITaxon:3983", "Manihot esculenta"),
    "chickpea": ("NCBITaxon:3827", "Cicer arietinum"),
    "soybean": ("NCBITaxon:3847", "Glycine max"),
    "soy": ("NCBITaxon:3847", "Glycine max"),
    "natane": ("NCBITaxon:3708", "Brassica napus"),
    "brassica napus": ("NCBITaxon:3708", "Brassica napus"),
    "allium": ("NCBITaxon:4678", "Allium"),
}

# Deliberately small: abbreviations, composite organs and organism-specific
# anatomical refinements require evidence beyond an apparent spelling match.
ORGANS = {
    "leaf": ("PO:0025034", "leaf"),
    "root": ("PO:0009005", "root"),
    "fruit": ("PO:0009001", "fruit"),
    "fruits": ("PO:0009001", "fruit"),
    "stem": ("PO:0009047", "stem"),
    "shoot": ("PO:0009006", "shoot system"),
    "whole plants": ("PO:0000003", "whole plant"),
    "whole plant": ("PO:0000003", "whole plant"),
    "seed": ("PO:0009010", "seed"),
    "seeds": ("PO:0009010", "seed"),
    "flower": ("PO:0009046", "flower"),
    "flowers": ("PO:0009046", "flower"),
    "ovary": ("PO:0009072", "plant ovary"),
    "pericarp": ("PO:0009084", "pericarp"),
    "nodule": ("PO:0003023", "root nodule"),
    "seedling": ("PO:0008037", "seedling"),
}

CHEMICALS = {
    "acoh": ("CHEBI:15366", "acetic acid"),
    "acetic acid": ("CHEBI:15366", "acetic acid"),
    "water": ("CHEBI:15377", "water"),
    "nacl": ("CHEBI:26710", "sodium chloride"),
    "peg": ("CHEBI:46793", "poly(ethylene glycol)"),
}

UNITS = {
    "mg": "UO:0000022", "g": "UO:0000021", "ng": "UO:0000024",
    "ml": "UO:0000098", "ul": "UO:0000101",
    "min": "UO:0000031", "h": "UO:0000032", "d": "UO:0000033",
    "w": "UO:0000034", "week": "UO:0000034", "weeks": "UO:0000034",
    "percent": "UO:0000187", "%": "UO:0000187",
    "um": "UO:0000064", "mm": "UO:0000063",
    "pmol": "UO:0000042", "cm2": "UO:0000081",
}
NUMBER = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"
QUANTITY_RE = re.compile(r"^(" + NUMBER + r")\s*([A-Za-z%µμ0-9^]*)$")


def _norm(value):
    return re.sub(r"\s+", " ", value.strip())


def _decimal(value):
    try:
        number = Decimal(value)
        if not number.is_finite():
            return None
        return format(number.normalize(), "f")
    except (InvalidOperation, ValueError):
        return None


@lru_cache(maxsize=1)
def _project_titles():
    """Use captured project titles only to disambiguate explicit name cells."""
    source = Path(__file__).resolve().parents[1] / "projects.json"
    try:
        return {item["p"].rsplit("/", 1)[-1]: re.sub(r"<[^>]*>", "", item.get("label", ""))
                for item in json.loads(source.read_text())}
    except (OSError, ValueError, KeyError):
        return {}


def _result(value, spec, **updates):
    result = {
        "value_kind": spec.get("value_kind", "unresolved"),
        "normalized_value": _norm(value),
        "term_id": "",
        "numeric_value": "",
        "unit_id": "",
        "reference_event": spec.get("reference_event", ""),
        "status": "needs_review",
        "reason": "No context-specific value mapping has been established.",
        "rule_id": "value.unresolved",
        "evidence": "source CSV column and value",
        "assertions": [],
    }
    result.update(updates)
    return result


def _term(value, spec, term_id, label, rule, evidence="source CSV column and value", **extra):
    fields = dict(value_kind="ontology_term", term_id=term_id,
                  normalized_value=label, status="accepted", rule_id=rule,
                  reason="Context-specific canonical term mapping; catalogue validation required.",
                  evidence=evidence)
    fields.update(extra)
    return _result(value, spec, **fields)


def _quantity(project, path, header, value, row, spec):
    val = _norm(value)
    low = val.lower()
    standard = spec.get("standard_name", "")
    unit = spec.get("unit_id", "")
    ref = spec.get("reference_event", "")
    reason = "Numeric literal with unit explicitly recorded in the column or cell."
    evidence = "source CSV column and value"
    status = "accepted"
    comparator = ""

    # Signed day offsets encode before/after flowering, not different units.
    if standard == "flowering_relative_time" or header.strip().lower() == "floweringday":
        match = re.fullmatch(r"(" + NUMBER + r")(dbf|daf)", low)
        if match:
            num = Decimal(match.group(1)) * (-1 if match.group(2) == "dbf" else 1)
            return _result(value, spec, value_kind="quantity", numeric_value=_decimal(str(num)),
                           unit_id="UO:0000033", reference_event="flowering", status="accepted",
                           rule_id="value.flowering_day_offset", reason="dbf/daf parsed as signed days before/after flowering.")

    if low.startswith("lessthan"):
        low = low[len("lessthan"):]
        comparator = "<"
        status = "needs_review"
        reason = "Recorded upper-bounded nominal mass class; not an exact measurement."

    dimensionless = standard in {"dilution_factor", "relative_dilution_factor", "dilution_ratio_description"}
    if dimensionless:
        stripped = re.sub(r"^con\s*:\s*", "", low)
        fraction = re.fullmatch(r"(" + NUMBER + r")\s*/\s*(" + NUMBER + r")", stripped)
        if fraction and Decimal(fraction.group(2)) != 0:
            number = Decimal(fraction.group(1)) / Decimal(fraction.group(2))
            return _result(value, spec, value_kind="quantity", numeric_value=_decimal(str(number)),
                           unit_id=unit, status="accepted", rule_id="value.dilution_fraction",
                           reason="Explicit dilution fraction parsed; original ratio retained in source value.")
        low = stripped

    match = QUANTITY_RE.fullmatch(low.replace("µ", "u").replace("μ", "u"))
    if not match:
        return _result(value, spec, rule_id="value.quantity_unparsed",
                       reason="Quantity cannot be parsed without interpreting an undocumented code.")
    num, suffix = match.groups()
    suffix_unit = UNITS.get(suffix)
    if suffix and suffix_unit:
        if unit and unit != suffix_unit:
            return _result(value, spec, value_kind="quantity", numeric_value=_decimal(num),
                           status="needs_review", rule_id="value.unit_conflict",
                           reason="The unit encoded in the value conflicts with the column unit.",
                           assertions=[{"role": "cell_unit", "term_id": suffix_unit},
                                       {"role": "column_unit", "term_id": unit}])
        unit = suffix_unit
    elif suffix:
        status = "needs_review"
        reason = "Numeric component parsed; the suffix has no verified unit mapping."

    # This file contains its own explicit mM NaCl evidence in the same row.
    if project == "RPMM0068" and header.strip().lower() == "conc":
        sample_name = next((v for k, v in row.items() if k.strip().lower() == "sample name"), "")
        explicit = re.search(r"\((" + NUMBER + r")\s*mM\s+NaCl\)", sample_name, re.I)
        if explicit and Decimal(explicit.group(1)) == Decimal(num):
            unit = "UO:0000063"
            evidence = "same-row sample name: " + sample_name
            reason = "Concentration unit and amount agree with explicit mM NaCl in the same-row sample name."
        else:
            status = "needs_review"
            reason = "The concentration needs a matching explicit unit from the same-row sample name."

    if not unit and not dimensionless:
        status = "needs_review"
        reason += " Unit is not explicitly established."
    if standard in {"elapsed_time", "recorded_elapsed_time", "flowering_relative_time", "flowering_related_period"} and not ref:
        status = "needs_review"
        reason += " Reference event for the elapsed time is not established."
    if spec.get("status") == "needs_review":
        status = "needs_review"
        reason += " Column semantics require review: " + spec.get("reason", "")
    assertions = []
    if comparator:
        assertions.append({"role": "numeric_comparator", "value": comparator})
    return _result(value, spec, value_kind="quantity", numeric_value=_decimal(num), unit_id=unit,
                   status=status, rule_id="value.quantity_literal", reason=reason,
                   evidence=evidence, assertions=assertions)


def annotate_value(project, path, header, value, row, spec):
    """Return one typed annotation plus optional structured assertions for a cell.

    ``row`` uses original CSV header names.  ``path`` is only used for explicit
    per-dataset rules, so both repository-relative and data-relative paths work.
    No result asserts that all rows in a project belong to one biological taxon.
    """
    value = "" if value is None else str(value)
    val = _norm(value)
    low = val.lower()
    role = spec.get("semantic_role", "")
    kind = spec.get("value_kind", "unresolved")
    standard = spec.get("standard_name", "")
    h = header.strip().lower()

    if not val:
        return _result(value, spec, value_kind="missing", status="missing", rule_id="value.empty",
                       reason="Empty source cell.")
    if low in {"nodata", "noinfo", "na", "n/a", "unknown", "unkown"}:
        return _result(value, spec, value_kind="missing", status="missing", rule_id="value.explicit_missing",
                       reason="Explicit unavailable/unknown marker; original spelling retained.")
    if low == "null":
        return _result(value, spec, value_kind="identifier", status="needs_review",
                       rule_id="value.null_genotype_retained",
                       reason="The genotype token null may denote a null segregant/no-insertion line; it is not converted to missing.")

    # Reference markers are not values for biological or physical properties.
    # FD still describes a biological reference material, so its missing
    # biological properties must not be declared inapplicable.  Notes and
    # protocol descriptions, however, genuinely contain these literal codes.
    if low in {"fd", "stds", "standard", "standards mixture", "standard mixture", "freeze dry col-0"}:
        if role in {"identifier", "file_reference", "order", "group"}:
            return _result(value, spec, value_kind="identifier", status="local", rule_id="value.local_reference_code",
                           reason="Dataset-local standard/reference code; not a universal biological class.")
        if role in {"text", "protocol_parameter", "compound_description"} or kind == "text":
            return _result(value, spec, value_kind="text", status="local", rule_id="value.reference_description",
                           reason="Reference/standard description or protocol code is retained as recorded text.")
        if low in {"stds", "standard", "standards mixture", "standard mixture"} and role in {"organism", "organism_part"}:
            return _result(value, spec, value_kind="text", status="not_applicable",
                           rule_id="value.standard_biological_property_inapplicable",
                           reason="Standard mixture marker, rather than a biological taxon or anatomical structure, occupies this source cell.")
        return _result(value, spec, value_kind="text", status="needs_review",
                       rule_id="value.reference_property_unknown",
                       reason="Reference/standard code does not state this biological property, quantity or treatment; applicability and value remain unestablished.")

    if kind == "quantity":
        return _quantity(project, path, header, value, row, spec)

    if kind == "integer" or role in {"order", "count"}:
        number = _decimal(val)
        if number is not None and Decimal(number) == Decimal(number).to_integral_value():
            return _result(value, spec, value_kind="integer", normalized_value=str(int(Decimal(number))),
                           numeric_value=number, status="accepted", rule_id="value.integer_literal",
                           reason="Recorded integer count or ordinal; the integer is not an ontology class.")
        return _result(value, spec, rule_id="value.integer_unparsed", reason="Expected an integer literal.")

    if kind == "boolean":
        if low in {"true", "false", "yes", "no"}:
            return _result(value, spec, value_kind="boolean", normalized_value="true" if low in {"true", "yes"} else "false",
                           status="accepted", rule_id="value.boolean_literal", reason="Explicit boolean state in a boolean column.")
        return _result(value, spec, rule_id="value.boolean_unparsed", reason="No documented boolean interpretation.")

    if kind == "date":
        if re.fullmatch(r"\d{8}", val):
            try:
                parsed = datetime.strptime(val, "%Y%m%d").date().isoformat()
                return _result(value, spec, value_kind="date", normalized_value=parsed, status="accepted",
                               rule_id="value.date_yyyymmdd", reason="Explicit eight-digit calendar date validated and normalized to ISO 8601.")
            except ValueError:
                pass
        if h == "samplingyear" and re.fullmatch(r"\d{4}", val):
            return _result(value, spec, value_kind="date", status="accepted", rule_id="value.calendar_year",
                           reason="Calendar year retained at year precision; no month/day is fabricated.")
        return _result(value, spec, value_kind="date", rule_id="value.ambiguous_date",
                       reason="Date format/century is not explicit enough to normalize without review.")

    if role == "organism":
        title = _project_titles().get(project, "")
        if low == "kanabikisou" and project == "RPMM0026" and "thesium chinense" in title.lower():
            return _term(value, spec, "NCBITaxon:210367", "Thesium chinense", "value.thesium_project_title",
                         evidence="explicit species=kanabikisou; investigation/projects.json " + project + ": " + title,
                         reason="Project title disambiguates the explicitly recorded common name; no taxon is inherited into other rows/columns.")
        if low == "arabidopsis":
            col0_fields = [(k, v) for k, v in row.items() if re.search(r"geno|line|tube|note", k, re.I)
                          and re.search(r"\bcol[-_ ]?0\b", str(v), re.I)]
            if "arabidopsis thaliana" in title.lower() or col0_fields:
                evidence = ("investigation/projects.json " + project + ": " + title
                            if "arabidopsis thaliana" in title.lower()
                            else "same-row Col-0 reference: " + "; ".join(k + "=" + str(v) for k, v in col0_fields))
                return _term(value, spec, "NCBITaxon:3702", "Arabidopsis thaliana", "value.arabidopsis_disambiguated",
                             evidence=evidence, reason="Explicit Arabidopsis name disambiguated to thaliana by the source title or same-row Col-0 reference.")
            return _result(value, spec, value_kind="ontology_term", rule_id="value.arabidopsis_genus_ambiguous",
                           reason="Arabidopsis is a genus name; no same-row Col-0 or explicit thaliana project-title evidence establishes the species.",
                           assertions=[{"role": "candidate_species", "term_id": "NCBITaxon:3702"}])
        if low == "natane" and "brassica napus" in title.lower():
            return _term(value, spec, "NCBITaxon:3708", "Brassica napus", "value.rapeseed_project_title",
                         evidence="explicit species=natane; investigation/projects.json " + project + ": " + title)
        if project == "RPMM0004" and h in {"speciesname", "binomialname"}:
            genus = _norm(row.get("material", "")).lower()
            epithet = _norm(row.get("binomialname", "")).lower()
            allium = {"fistulosum": ("NCBITaxon:35875", "Allium fistulosum"),
                      "chinense": ("NCBITaxon:130426", "Allium chinense"),
                      "tuberosum": ("NCBITaxon:4683", "Allium tuberosum")}
            if genus == "allium" and epithet in allium:
                term, label = allium[epithet]
                return _term(value, spec, term, label, "value.allium_genus_epithet",
                             evidence="same-row material=" + row.get("material", "") + "; binomialname=" + row.get("binomialname", ""))
        if standard == "substrate_or_inoculant_taxon" or project == "RPMM0037":
            # The taxon role changes from plant substrate to inoculated bacterium.
            # Keep the bacterial strain alongside a candidate taxon rather than
            # claiming that all measurements were made on a pure bacterial sample.
            if low == "soybean":
                return _term(value, spec, "NCBITaxon:3847", "Glycine max", "value.fermentation_substrate",
                             assertions=[{"role": "taxon_role", "value": "plant_substrate"}])
            microbial = re.fullmatch(r"bacillus (licheniformis|amyloliquefaciens)\s+(\S+)", low)
            if microbial:
                term = {"licheniformis": "NCBITaxon:1402", "amyloliquefaciens": "NCBITaxon:1390"}[microbial.group(1)]
                return _result(value, spec, value_kind="ontology_term", status="needs_review", rule_id="value.inoculant_taxon_strain",
                               reason="Named bacterial species and strain occur in a fermentation column; verify inoculant role and historical strain taxonomy.",
                               assertions=[{"role": "inoculant_taxon", "term_id": term},
                                           {"role": "inoculant_strain", "value": val.split()[-1]}])
        if low in SPECIES:
            term, label = SPECIES[low]
            return _term(value, spec, term, label, "value.explicit_organism_name")
        return _result(value, spec, value_kind="ontology_term", rule_id="value.ambiguous_organism_name",
                       reason="The common name/abbreviation does not establish a unique taxon without source metadata.")

    if role == "organism_part":
        if project == "RPMM0020" and low == "septum" and row.get("Plant", "").lower() == "tomato":
            return _term(value, spec, "PO:0000030", "septum", "value.tomato_fruit_septum",
                         evidence="same-row Plant=Tomato and Tissue=septum; investigation/projects.json RPMM0020 title identifies tomato fruit",
                         reason="Explicit septum in tomato-fruit tissue context matches the plant anatomical partition definition.")
        if low in ORGANS:
            term, label = ORGANS[low]
            if low == "seedling":
                return _term(value, spec, term, label, "value.explicit_seedling_material",
                             reason="Source seedling matches the native anatomical whole-plant class at the seedling development stage, not a developmental-stage class.",
                             assertions=[{"role": "sampled_material_kind", "value": "whole_plant_at_seedling_stage"}])
            return _term(value, spec, term, label, "value.explicit_plant_organ")
        if project == "RPMM0067" and low == "cortex" and row.get("tissue1", "").lower() == "root":
            return _term(value, spec, "PO:0000258", "root cortex", "value.root_cortex_hierarchy",
                         evidence="same-row tissue1=root; tissue3=Cortex")
        if project == "RPMM0067" and low == "parencyma":
            return _term(value, spec, "PO:0005421", "parenchyma", "value.parenchyma_curated_spelling",
                         evidence="RPMM0067 tissue3=Parencyma, inside root tissue hierarchy; curated spelling correction")
        if low in {"lb", "ls", "4lb"}:
            proposed = "PO:0020104" if low == "ls" else "PO:0020039"
            assertions = [{"role": "candidate_anatomical_structure", "term_id": proposed}]
            if low == "4lb":
                assertions.append({"role": "candidate_leaf_ordinal", "value": "4"})
            return _result(value, spec, value_kind="ontology_term", rule_id="value.anatomical_abbreviation_candidate",
                           reason="Candidate leaf lamina/sheath interpretation of source abbreviation needs confirmation from protocol/RDF.",
                           assertions=assertions)
        if low in {"septumseed", "jelly and seeds", "periexocarp", "noduleafe", "basalplate"}:
            return _result(value, spec, value_kind="ontology_term", rule_id="value.composite_or_abbreviated_organ",
                           reason="Composite, abbreviated or stage-specific anatomical value needs a verified decomposition/context.")
        return _result(value, spec, value_kind="ontology_term", rule_id="value.unmapped_plant_organ",
                       reason="No curated anatomical mapping exists for this spelling and tissue hierarchy.")

    if role == "genotype":
        if re.fullmatch(r"col[-_ ]?0", low):
            return _result(value, spec, value_kind="identifier", normalized_value="Col-0", status="local",
                           rule_id="value.col0_alias", reason="Col-0 spelling normalized; no stock accession is inferred.")
        return _result(value, spec, value_kind="identifier", status="local", rule_id="value.local_genotype",
                       reason="Study-local line/genotype identifier; an external accession has not been verified.")

    if role == "treatment":
        if standard == "growth_condition" and low == "soil":
            return _term(value, spec, "PECO:0007049", "soil environment exposure", "value.soil_growth_condition",
                         reason="Explicit soil value in growthcondition column identifies the growth medium; soil composition is not inferred.")
        if standard == "growth_condition" and low == "hc":
            return _result(value, spec, value_kind="ontology_term", rule_id="value.hydroponic_abbreviation_candidate",
                           reason="hc may denote hydroponic culture, but the abbreviation needs a protocol definition.",
                           assertions=[{"role": "candidate_growth_condition", "term_id": "PECO:0007067"}])
        if standard == "phosphorus_treatment" and ((project == "RPMM0044" and low == "p_stress") or
                                                    (project == "RPMM0062" and low == "minusp")):
            return _term(value, spec, "PECO:0001046", "limited phosphate exposure", "value.phosphate_deficiency_context",
                         evidence="source Ptreatment=" + val + "; investigation/projects.json " + project + " title explicitly describes phosphate deficiency/stress",
                         reason="Phosphorus treatment code denotes limitation in a project explicitly studying phosphate deficiency/stress; no dose is inferred.")
        if standard == "nitrogen_treatment" and low == "minusn":
            return _term(value, spec, "PECO:0001050", "limited nitrogen exposure", "value.explicit_nitrogen_withheld",
                         reason="Explicit minusN in Ntreatment records nitrogen withholding; its chemical form and dose are not inferred.")
        if low in CHEMICALS:
            term, label = CHEMICALS[low]
            return _term(value, spec, term, label, "value.explicit_treatment_chemical",
                         assertions=[{"role": "annotation_scope", "value": "chemical_agent_only"}],
                         reason="Chemical agent named in treatment cell; exposure protocol/condition is a separate annotation."
                                + (" Generic polymer only; molecular-weight grade/product is not inferred." if low == "peg" else ""))
        if low == "dry" and project == "RPMM0027":
            return _term(value, spec, "PECO:0007404", "drought exposure", "value.drought_condition",
                         evidence="source path RPMM0027/AnalysedDataset/DroughtTreatment; treatment=dry",
                         reason="Explicit dry code occurs in the DroughtTreatment dataset; records exposure without inferring severity or duration.")
        if low == "dex" and project == "RPMM0045":
            return _result(value, spec, value_kind="ontology_term", rule_id="value.dex_abbreviation_candidate",
                           reason="DEX commonly abbreviates dexamethasone in induction experiments; require the protocol to confirm reagent identity.",
                           assertions=[{"role": "candidate_chemical_agent", "term_id": "CHEBI:41879"}])
        if low in {"contrl", "control", "mock", "notreatment", "notreat", "nostress"}:
            return _result(value, spec, value_kind="text", normalized_value="control" if low == "contrl" else val,
                           status="local", rule_id="value.local_control_condition",
                           reason="Control/no-treatment code retained; its comparison and vehicle are study-specific.")
        return _result(value, spec, value_kind="ontology_term", rule_id="value.undocumented_treatment_code",
                       reason="Treatment code requires protocol context; chemical, stress and inoculation meanings are not interchangeable.")

    if role == "phenotype":
        colors = {"red": ("PATO:0000322", "red"), "pink": ("PATO:0000954", "pink")}
        if h == "color" and low in colors:
            term, label = colors[low]
            return _term(value, spec, term, label, "value.explicit_color")
        if h == "phenotype" and low == "elongated":
            return _term(value, spec, "PATO:0001154", "elongated", "value.explicit_shape",
                         reason="Explicit shape quality; the affected entity remains linked to the source row rather than inferred from the sampled tissue.")
        if h == "phenotype" and low == "round":
            return _result(value, spec, value_kind="ontology_term", rule_id="value.round_related_synonym",
                           reason="round is only a RELATED synonym of circular; dimensional/shape equivalence must be reviewed.",
                           assertions=[{"role": "candidate_shape_quality", "term_id": "PATO:0000411", "mapping_relation": "related"}])
        return _result(value, spec, value_kind="ontology_term", rule_id="value.phenotype_requires_context",
                       reason="Phenotype description needs an entity-quality or controlled trait mapping.")

    if role == "compound_description":
        return _result(value, spec, value_kind="text", status="needs_review", rule_id="value.compound_description",
                       reason="Multiple biological/experimental attributes are encoded together; source text retained for contextual decomposition.")

    if kind == "identifier" or role in {"identifier", "file_reference", "group", "inoculant"}:
        if h == "gene loci" and re.fullmatch(r"AT[1-5MC]G\d{5}", val, re.I):
            return _result(value, spec, value_kind="identifier", normalized_value=val.upper(), status="local",
                           rule_id="value.arabidopsis_locus", reason="Arabidopsis locus syntax recognized; no external database record/NCBI Gene accession is asserted.",
                           assertions=[{"role": "identifier_namespace", "value": "Arabidopsis locus nomenclature"}])
        return _result(value, spec, value_kind="identifier", status="local", rule_id="value.local_identifier",
                       reason="Local identifier, filename, group or strain code retained as a literal.")

    if kind == "text":
        return _result(value, spec, value_kind="text", status="local", rule_id="value.retained_text",
                       reason="Descriptive/protocol text retained without an unsupported ontology equivalence.")
    return _result(value, spec)
