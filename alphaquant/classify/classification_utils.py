
def find_node_parent_at_level(node, level):
    if node.level == level:
        return node
    while node.parent is not None:
        node = node.parent
        if node.level == level:
            return node

def find_node_parent_at_type(node, nodetype):
    if node.type == nodetype:
        return node
    while node.parent is not None:
        node = node.parent
        if node.type == nodetype:
            return node
