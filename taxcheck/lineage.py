import re

from functools import lru_cache

from ete3 import NCBITaxa


class Lineage:
	TAXLEVELS = {
        "kingdom": (0, "k"),
        "phylum": (1, "p"),
        "class": (2, "c"),
        "order": (3, "o"),
        "family": (4, "f"),
        "genus": (5, "g"),
        "species": (6, "s"),
    }

	def __init__(self, t_kingdom=None, t_phylum=None, t_class=None, t_order=None, t_family=None, t_genus=None, t_species=None):
		self.levels = (t_kingdom, t_phylum, t_class, t_order, t_family, t_genus, t_species)
	
	def get_string(self, show_names=True):
		key = 'name' if show_names else 'id'
		return ";".join(
			f"{prefix}__{self.levels[level][key]}"
			for level, prefix in Lineage.TAXLEVELS
			if self.levels[level] is not None
		)

	def get_lca(self, other):
		for i, (level_a, level_b) in enumerate(reversed(zip(self.levels, other.levels)), start=1):
			if level_a["id"] == level_b["id"]:
				return (7 - i), level_a
			
		return None


class LineageFactory:
	def __init__(self):
		self.ncbi = NCBITaxa()
	def generate_lineage(self, taxid):
		lineage = self.ncbi.get_lineage(taxid)
		ranks = {
			v:k
			for k,v in self.ncbi.get_rank(lineage).items()
		}
		ranks.setdefault("kingdom", ranks.get("superkingdom", None))
		names = self.ncbi.get_taxid_translator(lineage)
		return Lineage(**{
			rank: { "id": _id, "name": names.get(_id) } 
			for rank, _id in ranks.items()
			if Lineage.TAXLEVELS.get(rank)
		})



class LineageLookup:

    def __init__(self):
        self.ncbi = NCBITaxa()

    @lru_cache(maxsize=1000)
    def get_lineage(self, taxid):

        def get_lineage_str(levels):
            return ";".join(
                f"{prefix}__{string}"
                for _, prefix, string in sorted(levels)
            )

        try:
            lineage = self.ncbi.get_lineage(taxid)
        except ValueError:
            return (taxid, "INVALID", "INVALID")

        ranks = self.ncbi.get_rank(lineage)
        
        levels = []
        for tid, tname in zip(lineage, self.ncbi.translate_to_names(lineage)):
            rank = ranks.get(tid, "").replace("superkingdom", "kingdom")
            try:
                index, prefix = Lineage.TAXLEVELS.get(rank)
            except TypeError:
                continue
            levels.append((index, prefix, str(tid), re.sub(r" +", "_", tname)))

        return levels[-1][3], get_lineage_str((index, prefix, tid) for index, prefix, tid, _ in levels), get_lineage_str((index, prefix, tname) for index, prefix, _, tname in levels)




#>>> ncbi.get_lineage(662479)
#[1, 131567, 2157, 28890, 2290931, 183963, 1644055, 1644056, 2251, 403181, 662479]
#>>> ncbi.get_rank(ncbi.get_lineage(662479))
#{1: 'no rank', 2157: 'superkingdom', 2251: 'genus', 28890: 'phylum', 131567: 'no rank', 183963: 'class', 403181: 'species', 662479: 'strain', 1644055: 'order', 1644056: 'family', 2290931: 'clade'}
#>>> ncbi.translate_to_names(ncbi.get_lineage(662479))
#['root', 'cellular organisms', 'Archaea', 'Euryarchaeota', 'Stenosarchaea group', 'Halobacteria', 'Haloferacales', 'Haloferacaceae', 'Haloferax', 'Haloferax mucosum', 'Haloferax mucosum ATCC BAA-1512']
#>>> ncbi.get_taxid_translator(ncbi.get_lineage(662479))
#{1: 'root', 2157: 'Archaea', 2251: 'Haloferax', 28890: 'Euryarchaeota', 131567: 'cellular organisms', 183963: 'Halobacteria', 403181: 'Haloferax mucosum', 662479: 'Haloferax mucosum ATCC BAA-1512', 1644055: 'Haloferacales', 1644056: 'Haloferacaceae', 2290931: 'Stenosarchaea group'}
#
