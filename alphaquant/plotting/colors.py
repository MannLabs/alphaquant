import matplotlib

class AlphaQuantColorMap():
    def __init__(self):
        self.colorlist_hex  = ["#d8674e",  # Cadmium Red
        "#45a6ce",  # Steel Blue
        "#fdb73b",  # Cadmium Yellow
        "#a6d1f1",  # Baby Blue
        "#b04e8d",  # Tiffany Rose
        "#6e79b9",  # Periwinkle
        "#fcdf3b",  # Goldenrod
        "#50C878",  # Emerald Green
        "#808080",  # Grey instead of Amber
        "#FF7F50",  # Coral
        "#0F52BA",  # Egyptian Blue
        "#9966CC",  # Amethyst
        "#40E0D0"   # Turquoise
        ]
        self.colorlist = [matplotlib.colors.to_rgba(x) for x in self.colorlist_hex]


def rgba_list_to_hex_list(rgba_list):
    hex_list = []
    for rgba in rgba_list:
        # Convert each value to a 0-255 scale, then to hex, and finally concatenate.
        hex_code = '#' + ''.join([f"{int(c*255):02X}" for c in rgba[:3]])
        hex_list.append(hex_code)
    return hex_list


class AlphaPeptColorMap():
    def __init__(self):
        self.colorlist_hex  = ["#3FC5F0", "#42DEE1", "#7BEDC5", "#FFD479", "#16212B"]
        self.colorlist = [matplotlib.colors.to_rgba(x) for x in self.colorlist_hex]
        
        self.colormap_linear = matplotlib.colors.LinearSegmentedColormap.from_list("alphapept",self.colorlist)
        self.colormap_discrete = matplotlib.colors.LinearSegmentedColormap.from_list("alphapept",self.colorlist, N=5)

class ClusterColorMap():
    def __init__(self):
        self.colorlist_hex = [    "#D32F2F",  # Crimson Red
                    "#FFA000",  # Burnt Orange
                    "#FFEB3B",  # Golden Yellow
                    "#4CAF50",  # Grass Green
                    "#00BCD4",  # Cyan Blue
                    "#303F9F",  # Cobalt Blue
                    "#7B1FA2",  # Deep Purple
                    "#E91E63",  # Rose Pink
                    "#795548",  # Mocha Brown
                    "#607D8B"   # Slate Grey 
                    ]
        self.colorlist = [matplotlib.colors.to_rgba(x) for x in self.colorlist_hex]

class AlphaPeptColorMapAdapted():
    def __init__(self):
        self.colorlist_hex  =  [
    "#3FB7E4",  # Vivid Sky Blue (Slightly desaturated)
    "#7BEDC5",  # Medium Aquamarine
    "#EBCB70",  # Mustard (Slightly less bright)
    "#16212B",  # Gunmetal
    "#CA9ECB",  # Soft Lilac (Slightly more saturated)
    "#708090",  # Slate Gray
    "#3391A6",  # Deep Cerulean (Slightly lightened)
    "#AEDDE9",  # Powder Blue (Slightly warmer)
    "#5F9EA0",  # Cadet Blue
    "#E77D7D"   # Light Coral (Slightly desaturated)
]
        self.colorlist = [matplotlib.colors.to_rgba(x) for x in self.colorlist_hex]
        
