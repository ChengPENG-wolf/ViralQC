from pyvis.network import Network
import networkx as nx
import warnings


class Genome:

    def __init__(
        self,
        id: str,
        protein_info: list,
        protein_set: set,
        duplication_counts: dict,
    ):
        self.id = id
        self.protein_info = protein_info
        self.protein_set = protein_set
        self.duplication_counts = duplication_counts
        self.build_graph()

    def __add_gene(self, cluster, protein_id):
        self.graph.add_node(
            cluster,
            protein_id=protein_id,
        )

    def __add_edge(self, gene1, gene2):
        self.graph.add_edge(
            gene1,
            gene2,
            color='black',
            label=self.id
        )

    def build_graph(self):
        self.graph = nx.Graph(name=self.id)
        if len(self.protein_info) > 1:
            for i in range(0, len(self.protein_info) - 1):
                self.__add_gene(
                    self.protein_info[i][1],
                    self.protein_info[i][0],
                )
                self.__add_edge(
                    self.protein_info[i][1],
                    self.protein_info[i+1][1],
                )
            self.__add_gene(
                self.protein_info[-1][1],
                self.protein_info[-1][0],
            )

        elif len(self.protein_info) == 1:
            self.__add_gene(
                self.protein_info[0][1],
                self.protein_info[0][0],
            )
        else:
            warnings.warn(f'This genome ({self.id}) has no CDS.')

    def plot(self):
        nt = Network(width='1600px', height='1600px')
        nt.from_nx(self.graph)
        nt.show(f'{self.id}.html')


class GenomeStructure:
    def __init__(
        self,
        query_id,
        query_protein_info,
        ref_id,
        ref_protein_info,
        common_proteins,
    ):
        self.common_proteins = common_proteins
        self.events = {
            'match': 0,
            'mutation': 0,
            'insertion': 0,
            'deletion': 0,
            'translocation': 0,
            'duplication': 0,
            'outlier': 0,
            'gapopen': 0,
        }

        self.preprocess(query_id, query_protein_info, ref_id, ref_protein_info)
        self.get_match()
        self.get_indel()
        self.get_translocation()
        self.get_duplication()

    def preprocess(self, query_id, query_protein_info, ref_id, ref_protein_info):
        protein_set1 = set()
        dereplicated_protein_info1 = []
        duplication_counts1 = {}

        for protein in query_protein_info:
            cluster = protein[1]
            if cluster in protein_set1:
                duplication_counts1[cluster] = duplication_counts1.get(cluster, 0) + 1
            else:
                protein_set1.add(cluster)
                dereplicated_protein_info1.append(protein)

        protein_set2 = set()
        dereplicated_protein_info2 = []
        duplication_counts2 = {}

        for protein in ref_protein_info:
            cluster = protein[1]
            if cluster in protein_set2:
                duplication_counts2[cluster] = duplication_counts2.get(cluster, 0) + 1
            else:
                protein_set2.add(cluster)
                dereplicated_protein_info2.append(protein)

        self.unique_prots1 = protein_set1 - protein_set2
        self.unique_prots2 = protein_set2 - protein_set1
        self.unique_prots = self.unique_prots1 | self.unique_prots2

        # find the first and the last common proteins in dereplicated_protein_info1
        first = -1
        last = -1

        for i, protein in enumerate(dereplicated_protein_info1):
            if protein[1] in self.common_proteins:
                last = i
                if first == -1:
                    first = i
                    
        dereplicated_protein_info1 = dereplicated_protein_info1[first:last+1]
        self.genome1 = Genome(query_id, dereplicated_protein_info1, protein_set1, duplication_counts1)
        self.events['outlier'] += len(query_protein_info) - len(dereplicated_protein_info1)

        first = None
        last = None
        for i, protein in enumerate(dereplicated_protein_info2):
            if protein[1] in self.common_proteins:
                last = i
                if not first:
                    first = i
        dereplicated_protein_info2 = dereplicated_protein_info2[first:last + 1]
        self.genome2 = Genome(ref_id, dereplicated_protein_info2, protein_set2, duplication_counts2)

    def get_match(self):
        self.events['match'] = len(self.common_proteins)

    def get_indel(self):
        if len(self.common_proteins) == len(self.genome1.protein_set) == len(self.genome2.protein_set):
            return
        components = {}
        for genome in [self.genome1, self.genome2]:
            block = [genome.protein_info[0][1], None, []]

            for j in range(1, len(genome.protein_info) - 1):
                # if the protein is common, then the protein is the end of the block
                if genome.protein_info[j][1] in self.common_proteins:
                    # if the block is not empty, then add it to the components and reset the block
                    if len(block[2]) > 0:
                        block[1] = genome.protein_info[j][1]
                        if (block[0], block[1]) in components:
                            components[(block[0], block[1])].append([genome.id, block[2]])
                        elif (block[1], block[0]) in components:
                            components[(block[1], block[0])].append([genome.id, block[2]])
                        else:
                            components[(block[0], block[1])] = [[genome.id, block[2]]]

                        block = [None, None, []]

                    # if the next protein is not common, then the protein is the start of the block
                    if genome.protein_info[j + 1][1] not in self.common_proteins:
                        block[0] = genome.protein_info[j][1]
                # if the protein is not common, then add it to the block
                else:
                    block[2].append(genome.protein_info[j][1])

            # examine the final block
            if len(block[2]) > 0:
                block[1] = genome.protein_info[-1][1]
                if (block[0], block[1]) in components:
                    components[(block[0], block[1])].append([genome.id, block[2]])
                elif (block[1], block[0]) in components:
                    components[(block[1], block[0])].append([genome.id, block[2]])
                else:
                    components[(block[0], block[1])] = [[genome.id, block[2]]]

        # count the number of mutaitons, insertions, and deletions in each component
        for heads in components:
            component = components[heads]
            protein_g1 = []
            protein_g2 = []

            if len(component) == 1:
                if component[0][0] == self.genome1.id:
                    protein_g1 = component[0][1]
                else:
                    protein_g2 = component[0][1]
            else:
                protein_g1 = component[0][1]
                protein_g2 = component[1][1]

            self.events['mutation'] += min(len(protein_g1), len(protein_g2))
            if len(protein_g1) > len(protein_g2):
                self.events['insertion'] += (len(protein_g1) - len(protein_g2))
            else:
                self.events['deletion'] += (len(protein_g2) - len(protein_g1))

            if len(protein_g1) != len(protein_g2):
                self.events['gapopen'] += 1

    def get_translocation(self):
        if len(self.common_proteins) <= 2:
            return

        edges_g1 = set(self.genome1.graph.edges())
        edges_g2 = set(self.genome2.graph.edges())

        protein_info_g1 = [self.genome1.protein_info[0]]
        protein_info_g2 = [self.genome2.protein_info[0]]
        for i in range(1, len(self.genome1.protein_info) - 1):
            protein = self.genome1.protein_info[i]
            if protein[1] in self.common_proteins:
                protein_info_g1.append(protein)
            else:
                if (protein_info_g1[-1][1], protein[1]) in edges_g1:
                    edges_g1.remove((protein_info_g1[-1][1], protein[1]))
                else:
                    edges_g1.remove((protein[1], protein_info_g1[-1][1]))
                if (protein[1], self.genome1.protein_info[i+1][1]) in edges_g1:
                    edges_g1.remove((protein[1], self.genome1.protein_info[i+1][1]))
                else:
                    edges_g1.remove((self.genome1.protein_info[i+1][1], protein[1]))
                edges_g1.add((protein_info_g1[-1][1], self.genome1.protein_info[i+1][1]))

        for i in range(1, len(self.genome2.protein_info) - 1):
            protein = self.genome2.protein_info[i]
            if protein[1] in self.common_proteins:
                protein_info_g2.append(protein)
            else:
                if (protein_info_g2[-1][1], protein[1]) in edges_g2:
                    edges_g2.remove((protein_info_g2[-1][1], protein[1]))
                else:
                    edges_g2.remove((protein[1],  protein_info_g2[-1][1]))
                if (protein[1], self.genome2.protein_info[i+1][1]) in edges_g2:
                    edges_g2.remove((protein[1], self.genome2.protein_info[i+1][1]))
                else:
                    edges_g2.remove((self.genome2.protein_info[i+1][1], protein[1]))
                edges_g2.add((protein_info_g2[-1][1], self.genome2.protein_info[i+1][1]))

        self.tranloc_graph = nx.Graph(name=f'tranloc_{self.genome1.id}_{self.genome2.id}')

        for protein in self.common_proteins:
            self.tranloc_graph.add_node(
                protein
            )

        for a, b in edges_g1:
            if (a, b) in edges_g2 or (b, a) in edges_g2:
                self.tranloc_graph.add_edge(a, b, color='black')

        components = list(nx.connected_components(self.tranloc_graph))
        self.events['translocation'] = len(components) - 1

    def get_duplication(self):
        for protein in self.common_proteins:
            self.events['duplication'] += abs(self.genome1.duplication_counts.get(protein, 0) - self.genome2.duplication_counts.get(protein, 0))

    def get_events(self):
        return self.events

    def print_csv(self, filename):
        f = open(f'{filename}.csv', 'w')
        i = 0
        protein_id = {}
        for j, genome in enumerate([self.genome1, self.genome2]):
            for k in range(len(genome.protein_info) - 1):
                protein = genome.protein_info[k][1]
                if protein not in protein_id:
                    protein_id[protein] = i
                    i += 1
                if protein in self.common_proteins:
                    f.write(f'*{protein_id[protein]},')
                else:
                    f.write(f'{protein_id[protein]},')
            protein = genome.protein_info[-1][1]
            if protein not in protein_id:
                protein_id[protein] = i
                i += 1
            if protein in self.common_proteins:
                f.write(f'*{protein_id[protein]}\n')
            else:
                f.write(f'{protein_id[protein]}\n')
        f.write('\n')
        for protein in protein_id:
            f.write(f'{protein_id[protein]},{protein}\n')
        f.close()

    def get_joint_graph(self):
        self.joint_graph = nx.Graph(name=f'joint_{self.genome1.id}_{self.genome2.id}')

        for protein in self.genome1.protein_info:
            if protein[1] in self.common_proteins:
                self.joint_graph.add_node(
                    protein[1],
                    protein_id=protein[0],
                )
            else:
                self.joint_graph.add_node(
                    protein[1],
                    protein_id=protein[0],
                    shape='square',
                    color='#bdc0e8'
                )

        for protein in self.genome2.protein_info:
            if protein[1] not in self.common_proteins:
                self.joint_graph.add_node(
                    protein[1],
                    protein_id=protein[0],
                    shape='triangle',
                    color='#fbc4d6'
                )

        for a, b in self.genome1.graph.edges():
            if self.genome2.graph.has_edge(a, b):
                self.joint_graph.add_edge(a, b, color='black')
            else:
                self.joint_graph.add_edge(a, b, color='blue')

        for a, b in self.genome2.graph.edges():
            if not self.genome1.graph.has_edge(a, b):
                self.joint_graph.add_edge(a, b, color='red')

    def plot_joint_graph(self):
        nt = Network(width='1600px', height='1600px')
        nt.from_nx(self.joint_graph)
        nt.show(f'joint_{self.genome1.id}_{self.genome2.id}.html')

    def get_simplified_joint_graph(self):
        self.s_joint_graph = nx.Graph(name=f's_joint_{self.genome1.id}_{self.genome2.id}')
        self.s_joint_graph.add_nodes_from(self.joint_graph.nodes(data=True))
        self.s_joint_graph.add_edges_from(self.joint_graph.edges(data=True))

        for node in self.joint_graph.nodes():
            node_edges = [(a, b, att) for a, b, att in self.s_joint_graph.edges(data=True) if node in (a, b)]
            important = False
            for edge in node_edges:
                if edge[2]['color'] in ['blue', 'red']:
                    important = True
                    break
            if not important:
                node_neighbors = list(self.s_joint_graph.neighbors(node))
                self.s_joint_graph.remove_node(node)
                for i, neighbor1 in enumerate(node_neighbors):
                    for j, neighbor2 in enumerate(node_neighbors):
                        if i < j:
                            self.s_joint_graph.add_edge(neighbor1, neighbor2, color='black')

    def plot_simplified_joint_graph(self):
        nt = Network(width='1600px', height='1600px')
        nt.from_nx(self.s_joint_graph)
        nt.show(f's_joint_{self.genome1.id}_{self.genome2.id}.html')

    def plot_translocation(self):
        nt = Network(width='1600px', height='1600px')
        nt.from_nx(self.tranloc_graph)
        nt.show(f'tranloc_{self.genome1.id}_{self.genome2.id}.html')





