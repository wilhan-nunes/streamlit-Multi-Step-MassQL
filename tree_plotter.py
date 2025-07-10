import plotly.graph_objects as go
from utils import bile_acid_tree


class BileAcidTreeVisualizer:
    def __init__(self, tree_dict):
        self.tree_dict = tree_dict

    def build_sankey_data(self, highlight_path=None):
        """Build nodes and links for Sankey diagram with optional path highlighting."""
        nodes = []
        links = []
        node_dict = {}
        node_counter = 0

        def add_nodes_and_links(data, parent_idx=None, level=0):
            nonlocal node_counter

            for key, value in data.items():
                # Add node
                current_idx = node_counter
                node_dict[key] = current_idx

                # Color based on highlighting and level
                if highlight_path and key in highlight_path:
                    color = 'rgba(255, 107, 107, 0.9)'  # Bright red for highlighted path
                else:
                    # Different colors for different levels
                    level_colors = [
                        'rgba(78, 205, 196, 0.8)',   # Teal - root level
                        'rgba(255, 183, 77, 0.8)',   # Orange - stage level
                        'rgba(129, 236, 236, 0.8)',  # Light blue - intermediate
                        'rgba(255, 218, 185, 0.8)',  # Peach - specific compounds
                        'rgba(162, 155, 254, 0.8)',  # Purple - final level
                    ]
                    color = level_colors[min(level, len(level_colors)-1)]

                nodes.append({
                    'label': key,
                    'color': color
                })
                node_counter += 1

                # Add link from parent with highlighting
                if parent_idx is not None:
                    # Highlight links that are part of the path
                    if highlight_path and key in highlight_path:
                        parent_key = [k for k, v in node_dict.items() if v == parent_idx][0]
                        if parent_key in highlight_path:
                            link_color = 'rgba(255, 107, 107, 0.8)'  # Red for highlighted links
                        else:
                            link_color = 'rgba(200, 200, 200, 0.4)'  # Gray for non-highlighted
                    else:
                        link_color = 'rgba(200, 200, 200, 0.4)'  # Default gray

                    links.append({
                        'source': parent_idx,
                        'target': current_idx,
                        'value': 1,
                        'color': link_color
                    })

                # Process children recursively
                if isinstance(value, dict) and value:
                    add_nodes_and_links(value, current_idx, level + 1)

        add_nodes_and_links(bile_acid_tree)
        return nodes, links, node_dict

    def create_sankey_diagram(self, highlight_path=None, title_suffix=""):
        """Create a Sankey diagram with optional path highlighting."""
        nodes, links, node_dict = self.build_sankey_data(highlight_path)
        
        fig = go.Figure(data=[go.Sankey(
            node=dict(
                pad=15,
                thickness=20,
                line=dict(color="black", width=0.5),
                label=[node['label'] for node in nodes],
                color=[node['color'] for node in nodes],
                align="left",
            ),
            link=dict(
                arrowlen=15,
                source=[link['source'] for link in links],
                target=[link['target'] for link in links],
                value=[link['value'] for link in links],
                color=[link['color'] for link in links]
            ),
            textfont=dict(
                size=16,
                color="black",
                shadow="0px -0px 2px white"),
        )])
        
        title = f"Bile Acid Classification Tree{title_suffix}"
        if highlight_path:
            title += f"<br><sub>Highlighted Path: {' → '.join(highlight_path)}</sub>"
        
        fig.update_layout(
            title=title,
            font_size=14,
            width=1400,
            height=800,
            # title_x=0.5
        )
        
        return fig
    
    def create_multiple_diagrams(self, paths_dict):
        """Create multiple diagrams showing different highlighted paths."""
        figs = []
        
        # Create base diagram without highlighting
        base_fig = self.create_sankey_diagram(title_suffix=" - Base View")
        figs.append(('Base View', base_fig))
        
        # Create highlighted diagrams for each path
        for path_name, path in paths_dict.items():
            fig = self.create_sankey_diagram(path, f" - {path_name}")
            figs.append((path_name, fig))
        
        return figs

# Initialize the visualizer
visualizer = BileAcidTreeVisualizer(bile_acid_tree)

# Create specific diagram with custom path
def create_custom_tree(custom_path, path_name="Custom Path"):
    """Create diagram with custom highlighted path."""
    fig = visualizer.create_sankey_diagram(custom_path, f" - {path_name}")
    return fig


# Function to easily highlight any path
def highlight_path(path, name="Highlighted Path"):
    """Convenience function to highlight a specific path."""
    fig = visualizer.create_sankey_diagram(path, f" - {name}")
    fig.show()
    return fig


## Example usage:
# Print available paths for reference
# print("Available predefined paths:")
# for name, path in example_paths.items():
#     print(f"- {name}: {' → '.join(path)}")
