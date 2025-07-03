from collections import deque

class Edge:
    def __init__(self, source=None, target=None, weight=1):
        self.source = source
        self.target = target
        self.weight = weight

    def __str__(self):
        return f"({self.source.id} -> {self.target.id}, {self.weight})"

    def __repr__(self):
        return f"({self.source.id} -> {self.target.id}, {self.weight})"



class Node:
    def __init__(self, id, value=None):
        self.id = id
        self.value = value
        self.in_edges = []
        self.out_edges = []

    def __str__(self):
        return f"id: {self.id}; out edges: {[repr(edge) for edge in self.out_edges]}"

    def __repr__(self):
        return f"id: {self.id}; out edges: {[repr(edge) for edge in self.out_edges]}\n"



class Graph:
    def __init__(self, is_directed=True):
        self.nodes = {}
        self.is_directed = is_directed

    def add_node(self, id, value=None):
        if id not in self.nodes:
            self.nodes[id] = Node(id, value)

    def add_edge(self, source_id, target_id, weight=1):
        if (source_id not in self.nodes) or (target_id not in self.nodes):
            raise Exception('at least one node does not exist in graph')

        source = self.nodes[source_id]
        target = self.nodes[target_id]

        edge = Edge(source, target, weight)
        source.out_edges.append(edge)
        target.in_edges.append(edge)

    def get_node(self, id):
        return self.nodes[id]

    def __str__(self):
        return f"{[repr(node) for node in self.nodes.values()]}"



    def topological_ordering(self):
        sorted_nodes = []
        in_degrees = {}
        queue = deque()
        for node in self.nodes.values():
            in_degrees[node.id] = len(node.in_edges)
            if in_degrees[node.id] == 0:
                queue.append(node)

        while len(queue) > 0:
            current_node = queue.popleft()
            sorted_nodes.append(current_node)
            for edge in current_node.out_edges:
                target = edge.target
                in_degrees[target.id] -= 1
                if in_degrees[target.id] == 0:
                    queue.append(target)

        if len(sorted_nodes) < len(self.nodes):
            raise Exception('graph contains cycle')

        return sorted_nodes



    def longest_path_DAG(self, source_id, sink_id):
        scores = {} # stores tuples (score, backtracking_id)
        for node_id in self.nodes:
            scores[node_id] = (float('-inf'), None)
        scores[source_id] = (0, None)

        sorted_graph = self.topological_ordering()
        sorted_graph = sorted_graph[1:] # remove the source node
        for node in sorted_graph:
            max_score_tuple = (float('-inf'), None)
            for edge in node.in_edges:
                score = scores[edge.source.id][0] + edge.weight
                if score > max_score_tuple[0]:
                    max_score_tuple = (score, edge.source.id)
            scores[node.id] = max_score_tuple

        # backtracking
        longest_path = []
        temp_node_id = sink_id
        while temp_node_id != source_id:
            longest_path.append(temp_node_id)
            temp_node_id = scores[temp_node_id][1]
        longest_path.append(source_id)

        longest_path.reverse()
        return longest_path




