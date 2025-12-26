def validate_consistency(tree, strict=True):
    errs = []

    def err(path, msg):
        errs.append(f"{path}: {msg}")

    def walk(node, path, parent):
        s = node["start"]; e = node["end"]; L = node["length"]
        if s is not None or e is not None:
            if not isinstance(s,int) or not isinstance(e,int):
                err(path, "start/end must be int when present")
            else:
                if s > e: err(path, "start > end")
                if L is None or L != e - s + 1:
                    err(path, f"length mismatch (length={L}, expected={e - s + 1})")

        if L is not None and node["missing_N_count"] > L:
            err(path, "missing_N_count > length")

        if parent and parent["start"] is not None and parent["end"] is not None and s is not None and e is not None:
            if s < parent["start"] or e > parent["end"]:
                err(path, "child interval out of parent bounds")

        # gene partition checks
        if node["type"] == "gene":
            exon_g = None; intron_g = None
            for c in node["children"]:
                if c["type"] == "exon_group": exon_g = c
                if c["type"] == "non_exon_group": intron_g = c
            if exon_g and intron_g and node["length"] is not None:
                if (exon_g["length"] or 0) + (intron_g["length"] or 0) != node["length"]:
                    err(path, "exon+non_exon length != gene length")
                if exon_g["variant_count"] + intron_g["variant_count"] != node["variant_count"]:
                    err(path, "exon+non_exon variants != gene variants")

        # exon partition checks if exon_other_group exists
        if node["type"] == "exon_segment":
            cds_g = utr_g = oth_g = None
            for c in node["children"]:
                if c["type"] == "cds_group": cds_g = c
                if c["type"] == "utr_group": utr_g = c
                if c["type"] == "exon_other_group": oth_g = c
            if cds_g and utr_g and oth_g and node["length"] is not None:
                if (cds_g["length"] or 0) + (utr_g["length"] or 0) + (oth_g["length"] or 0) != node["length"]:
                    err(path, "cds+utr+other length != exon length")
                if cds_g["variant_count"] + utr_g["variant_count"] + oth_g["variant_count"] != node["variant_count"]:
                    err(path, "cds+utr+other variants != exon variants")

        for i, c in enumerate(node["children"]):
            walk(c, f"{path}.children[{i}]", node)

    walk(tree, "root", None)
    return errs
