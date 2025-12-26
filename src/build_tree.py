from dataclasses import dataclass
from pathlib import Path
import statistics

from utils.discovery import discover_assembly_files
from utils.intervals import merge_intervals, complement_intervals, intersect_intervals, subtract_intervals, interval_total_len, count_hit
from parsers.fasta_scan import scan_fasta_files
from parsers.gff_ensembl import parse_gff_files
from parsers.vcf_stream import iter_vcf_variants
from tree.node import RegionNode

@dataclass
class GeneRuntime:
    gene: RegionNode
    exon_segs: list  # [(s,e,node)]
    intron_segs: list
    exon_sub: dict   # exon_node -> {"cds":[(s,e,node)],"utr":[...],"other":[...]}

def _n_overlap(n_runs, interval):
    s,e = interval
    cnt = 0
    runs = 0
    i = 0
    while i < len(n_runs) and n_runs[i][1] < s:
        i += 1
    while i < len(n_runs) and n_runs[i][0] <= e:
        ns, ne = n_runs[i]
        is_ = max(s, ns); ie = min(e, ne)
        if is_ <= ie:
            cnt += (ie - is_ + 1)
            runs += 1
        i += 1
    return cnt, runs

def _make_group(node_type, base_id, chrom, segs, n_runs):
    group = RegionNode(node_type=node_type, node_id=base_id, chrom=chrom)
    seg_nodes = []
    for i,(s,e) in enumerate(segs, start=1):
        seg_type = node_type.replace("_group","_segment")
        seg = RegionNode(node_type=seg_type, node_id=f"{base_id}:{i}", chrom=chrom, start=s, end=e)
        ncnt, nruns = _n_overlap(n_runs, (s,e))
        seg.missing_N_count = ncnt
        seg.missing_N_runs = nruns
        group.add(seg)
        seg_nodes.append((s,e,seg))
    group.length = interval_total_len(segs)
    group.missing_N_count = sum(n.missing_N_count for _,_,n in seg_nodes)
    group.missing_N_runs = sum(n.missing_N_runs for _,_,n in seg_nodes)
    return group, seg_nodes

def build_species_tree(species_dir: Path, strict: bool=True) -> RegionNode:
    species = RegionNode(node_type="species", node_id=species_dir.name)

    for assembly_dir in sorted([p for p in species_dir.iterdir() if p.is_dir()]):
        files = discover_assembly_files(assembly_dir)
        asm = RegionNode(node_type="assembly", node_id=assembly_dir.name)
        asm.meta["file_counts"] = {"fasta": len(files.fasta), "gff": len(files.gff), "vcf": len(files.vcf)}
        species.add(asm)

        if not files.fasta:
            if strict:
                raise RuntimeError(f"[{species_dir.name}/{assembly_dir.name}] No FASTA files found")
            else:
                continue

        fasta = scan_fasta_files(files.fasta)
        gff = parse_gff_files(files.gff) if files.gff else {}

        gene_runtime = {}
        chrom_nodes = {}

        for chrom, finfo in fasta.items():
            chr_node = RegionNode(node_type="chromosome", node_id=chrom, chrom=chrom, start=1, end=finfo.length)
            chr_node.missing_N_count = finfo.missing_N_count
            chr_node.missing_N_runs = finfo.missing_N_runs
            asm.add(chr_node)
            chrom_nodes[chrom] = chr_node

            genes_on_chr = gff.get(chrom, {})
            gene_runtime[chrom] = {}

            gene_union = RegionNode(node_type="gene_union_group", node_id=f"{chrom}:gene_union", chrom=chrom, start=1, end=finfo.length)
            non_gene = RegionNode(node_type="non_gene_group", node_id=f"{chrom}:non_gene", chrom=chrom, start=1, end=finfo.length)
            chr_node.add(gene_union)
            chr_node.add(non_gene)

            gene_spans = merge_intervals([(gm.start, gm.end) for gm in genes_on_chr.values()])
            non_gene_segs = complement_intervals((1, finfo.length), gene_spans)
            ng_group, ng_nodes = _make_group("non_gene_group", f"{chrom}:non_gene_seg", chrom, non_gene_segs, finfo.n_runs)
            non_gene.children = ng_group.children
            non_gene.length = ng_group.length
            non_gene.missing_N_count = ng_group.missing_N_count
            non_gene.missing_N_runs = ng_group.missing_N_runs

            for gid, gm in genes_on_chr.items():
                gene = RegionNode(node_type="gene", node_id=gid, chrom=chrom, start=gm.start, end=gm.end)
                if gm.biotype:
                    gene.meta["biotype"] = gm.biotype
                ncnt, nruns = _n_overlap(finfo.n_runs, (gm.start, gm.end))
                gene.missing_N_count = ncnt
                gene.missing_N_runs = nruns
                gene_union.add(gene)

                exons = merge_intervals([(max(gm.start,s), min(gm.end,e)) for s,e in gm.exon if not (e < gm.start or s > gm.end)])
                introns = complement_intervals((gm.start, gm.end), exons)

                exon_group, exon_nodes = _make_group("exon_group", f"{gid}:exon", chrom, exons, finfo.n_runs)
                intron_group, intron_nodes = _make_group("non_exon_group", f"{gid}:intron", chrom, introns, finfo.n_runs)
                gene.add(exon_group)
                gene.add(intron_group)

                cds = merge_intervals([(max(gm.start,s), min(gm.end,e)) for s,e in gm.cds if not (e < gm.start or s > gm.end)])
                utr = merge_intervals([(max(gm.start,s), min(gm.end,e)) for s,e in gm.utr if not (e < gm.start or s > gm.end)])
                exon_sub = {}

                for s,e,enode in exon_nodes:
                    cds_in = merge_intervals(intersect_intervals(cds, [(s,e)]))
                    utr_in = merge_intervals(intersect_intervals(utr, [(s,e)]))
                    other_in = subtract_intervals([(s,e)], merge_intervals(cds_in + utr_in))

                    cds_g, cds_nodes = _make_group("cds_group", f"{enode.node_id}:cds", chrom, cds_in, finfo.n_runs)
                    utr_g, utr_nodes = _make_group("utr_group", f"{enode.node_id}:utr", chrom, utr_in, finfo.n_runs)
                    oth_g, oth_nodes = _make_group("exon_other_group", f"{enode.node_id}:other", chrom, other_in, finfo.n_runs)

                    enode.add(cds_g); enode.add(utr_g); enode.add(oth_g)

                    exon_sub[enode] = {"cds": cds_nodes, "utr": utr_nodes, "other": oth_nodes}

                gene_runtime[chrom][gid] = GeneRuntime(gene, exon_nodes, intron_nodes, exon_sub)

        # stream variants
        gene_spans_by_chr = {}
        for chrom, genes in gene_runtime.items():
            spans = [(rt.gene.start, rt.gene.end, gid) for gid, rt in genes.items()]
            spans.sort(key=lambda x: (x[0], x[1]))
            gene_spans_by_chr[chrom] = spans

        active_by_chr = {c: [] for c in gene_spans_by_chr}
        idx_by_chr = {c: 0 for c in gene_spans_by_chr}

        # non-gene segment lists
        non_gene_segs_by_chr = {}
        for chrom, chr_node in chrom_nodes.items():
            ng = None
            for c in chr_node.children:
                if c.node_type == "non_gene_group":
                    ng = c
                    break
            if ng:
                segs = [(ch.start, ch.end, ch) for ch in ng.children if ch.start and ch.end]
                segs.sort(key=lambda x: (x[0], x[1]))
                non_gene_segs_by_chr[chrom] = segs

        for chrom, pos in iter_vcf_variants(files.vcf):
            if chrom not in chrom_nodes:
                if strict:
                    raise RuntimeError(f"[{species_dir.name}/{assembly_dir.name}] VCF chrom not in FASTA: {chrom}")
                else:
                    continue
            chr_node = chrom_nodes[chrom]
            if pos < 1 or pos > chr_node.end:
                if strict:
                    raise RuntimeError(f"[{species_dir.name}/{assembly_dir.name}] VCF pos out of bounds: {chrom}:{pos} len={chr_node.end}")
                else:
                    continue

            chr_node.variant_count += 1

            # non-gene
            hit = count_hit(non_gene_segs_by_chr.get(chrom, []), pos)
            if hit is not None:
                hit.variant_count += 1

            spans = gene_spans_by_chr.get(chrom, [])
            k = idx_by_chr.get(chrom, 0)
            active = active_by_chr.get(chrom, [])

            while k < len(spans) and spans[k][0] <= pos:
                active.append(spans[k])  # (s,e,gid)
                k += 1
            idx_by_chr[chrom] = k

            active = [t for t in active if t[1] >= pos]
            active_by_chr[chrom] = active

            for gs, ge, gid in active:
                if gs <= pos <= ge:
                    rt = gene_runtime[chrom].get(gid)
                    if not rt:
                        continue
                    rt.gene.variant_count += 1

                    exon_hit = count_hit(rt.exon_segs, pos)
                    if exon_hit is not None:
                        exon_hit.variant_count += 1
                        sub = rt.exon_sub.get(exon_hit, {})
                        for key in ("cds","utr","other"):
                            node = count_hit(sub.get(key, []), pos)
                            if node is not None:
                                node.variant_count += 1
                    else:
                        node = count_hit(rt.intron_segs, pos)
                        if node is not None:
                            node.variant_count += 1

        # compute derived lengths
        species.compute_lengths()

        # aggregate group nodes' variant counts as sum(children)
        def aggregate(node: RegionNode):
            for c in node.children:
                aggregate(c)
            if node.node_type.endswith("_group") or node.node_type in ("species","assembly"):
                node.variant_count = sum(ch.variant_count for ch in node.children)
                if node.node_type not in ("chromosome","gene"):
                    node.missing_N_count = sum(ch.missing_N_count for ch in node.children)
                    node.missing_N_runs = sum(ch.missing_N_runs for ch in node.children)
                if node.length is None:
                    node.length = sum((ch.length or 0) for ch in node.children)
        aggregate(species)

        # assembly gene distributions
        gene_lengths = []
        gene_vars = []
        for chrom, genes in gene_runtime.items():
            for gid, rt in genes.items():
                if rt.gene.start is not None and rt.gene.end is not None:
                    gene_lengths.append(rt.gene.end - rt.gene.start + 1)
                    gene_vars.append(rt.gene.variant_count)
        if gene_lengths:
            asm.meta["gene_length_dist"] = {
                "count": len(gene_lengths),
                "min": min(gene_lengths),
                "max": max(gene_lengths),
                "mean": statistics.mean(gene_lengths),
                "var": statistics.pvariance(gene_lengths) if len(gene_lengths) > 1 else 0.0,
            }
            asm.meta["gene_variant_dist"] = {
                "count": len(gene_vars),
                "min": min(gene_vars),
                "max": max(gene_vars),
                "mean": statistics.mean(gene_vars),
                "var": statistics.pvariance(gene_vars) if len(gene_vars) > 1 else 0.0,
            }

    return species
