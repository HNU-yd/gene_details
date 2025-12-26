REQUIRED_KEYS = {"type","id","chrom","start","end","length","variant_count","missing_N_count","missing_N_runs","meta","children"}

def validate_structure(node, path="root"):
    errs = []
    if not isinstance(node, dict):
        return [f"{path}: not a dict"]
    missing = REQUIRED_KEYS - set(node.keys())
    if missing:
        errs.append(f"{path}: missing keys {sorted(missing)}")
        return errs

    if not isinstance(node["type"], str): errs.append(f"{path}.type not str")
    if not isinstance(node["id"], str): errs.append(f"{path}.id not str")
    if not (node["chrom"] is None or isinstance(node["chrom"], str)): errs.append(f"{path}.chrom invalid")
    if not (node["start"] is None or isinstance(node["start"], int)): errs.append(f"{path}.start invalid")
    if not (node["end"] is None or isinstance(node["end"], int)): errs.append(f"{path}.end invalid")
    if not (node["length"] is None or isinstance(node["length"], int)): errs.append(f"{path}.length invalid")
    for k in ("variant_count","missing_N_count","missing_N_runs"):
        if not isinstance(node[k], int) or node[k] < 0:
            errs.append(f"{path}.{k} invalid")
    if not isinstance(node["meta"], dict): errs.append(f"{path}.meta not dict")
    if not isinstance(node["children"], list): errs.append(f"{path}.children not list")

    for i, c in enumerate(node["children"]):
        errs.extend(validate_structure(c, f"{path}.children[{i}]"))
    return errs
