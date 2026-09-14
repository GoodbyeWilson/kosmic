"""
Pathway Registry
=================
Central in-memory store for gene set collections with search, coverage,
and import/export capabilities.
Pure Python — no GUI imports.
"""

from kosmic.reference.pathways import formats


class PathwayRegistry:
    """In-memory pathway store with collections, search, coverage, and I/O."""

    def __init__(self):
        self.collections = {}  # {collection_name: {pathway_name: [genes]}}

    # ------------------------------------------------------------------
    # Collection management
    # ------------------------------------------------------------------

    def add_collection(self, name, pathways):
        """Add or replace a named collection of pathways.

        Parameters
        ----------
        name : str
            Collection name (e.g. 'Metabolic (Comprehensive)', 'KEGG_2026').
        pathways : dict
            {pathway_name: [gene_list]}.
        """
        self.collections[name] = dict(pathways)

    def get_collection(self, name):
        """Return a copy of a single collection's pathways.

        Returns
        -------
        dict or None
            {pathway_name: [genes]} or None if not found.
        """
        coll = self.collections.get(name)
        return dict(coll) if coll is not None else None

    # ------------------------------------------------------------------
    # Merged / flat access
    # ------------------------------------------------------------------

    def get_flat(self):
        """Merge all collections into a single {pathway: [genes]} dict.

        If the same pathway name appears in multiple collections, the
        gene lists are merged (union).

        Returns
        -------
        dict
            {pathway_name: sorted [gene_list]}.
        """
        merged = {}
        for coll in self.collections.values():
            for pathway, genes in coll.items():
                if pathway in merged:
                    merged[pathway] = sorted(set(merged[pathway]) | set(genes))
                else:
                    merged[pathway] = list(genes)
        return merged

    def total_genes(self):
        """Count total unique genes across all collections."""
        genes = set()
        for coll in self.collections.values():
            for gene_list in coll.values():
                genes.update(gene_list)
        return len(genes)

    # ------------------------------------------------------------------
    # Search
    # ------------------------------------------------------------------

    def search(self, query):
        """Search pathway names across all collections (case-insensitive).

        Parameters
        ----------
        query : str
            Search term.

        Returns
        -------
        list of tuple
            [(collection_name, pathway_name, gene_count), ...] matching results.
        """
        query_lower = query.lower()
        results = []
        for coll_name, coll in self.collections.items():
            for pathway, genes in coll.items():
                if query_lower in pathway.lower():
                    results.append((coll_name, pathway, len(genes)))
        return results

    # ------------------------------------------------------------------
    # Coverage
    # ------------------------------------------------------------------

    # ------------------------------------------------------------------
    # File I/O
    # ------------------------------------------------------------------

    def load_file(self, path, collection_name=None):
        """Load gene sets from a file and add as a collection.

        Parameters
        ----------
        path : str or Path
            Path to .gmt, .json, or .csv file.
        collection_name : str, optional
            Name for the collection. Defaults to the filename stem.

        Returns
        -------
        str
            The collection name used.
        """
        from pathlib import Path as P
        p = P(path)
        name = collection_name or p.stem
        pathways = formats.load_file(p)
        self.add_collection(name, pathways)
        return name

    def load_directory(self, path):
        """Load every .gmt / .json / .csv file in a directory as a collection.

        The filename stem becomes the collection name. Subdirectories are
        ignored.

        Parameters
        ----------
        path : str or Path
            Directory containing gene set files.

        Returns
        -------
        list of str
            Collection names that were loaded.
        """
        from pathlib import Path as P
        directory = P(path)
        loaded = []
        for ext in ('*.gmt', '*.json', '*.csv'):
            for f in sorted(directory.glob(ext)):
                loaded.append(self.load_file(f))
        return loaded

    def save_gmt(self, path, collection_name=None):
        """Save pathways to GMT format.

        Parameters
        ----------
        path : str or Path
            Output path.
        collection_name : str, optional
            Save a specific collection. If None, saves the merged flat set.
        """
        if collection_name:
            pathways = self.collections.get(collection_name, {})
        else:
            pathways = self.get_flat()
        formats.save_gmt(pathways, path)

    def save_json(self, path, collection_name=None):
        """Save pathways to JSON format.

        Parameters
        ----------
        path : str or Path
            Output path.
        collection_name : str, optional
            Save a specific collection. If None, saves the merged flat set.
        """
        if collection_name:
            pathways = self.collections.get(collection_name, {})
        else:
            pathways = self.get_flat()
        formats.save_json(pathways, path)

    # ------------------------------------------------------------------
    # String representation
    # ------------------------------------------------------------------

    def __repr__(self):
        parts = []
        for name, coll in self.collections.items():
            n_genes = len({g for genes in coll.values() for g in genes})
            parts.append(f"  {name}: {len(coll)} pathways, {n_genes} genes")
        body = '\n'.join(parts) if parts else '  (empty)'
        return f"PathwayRegistry(\n{body}\n)"
