import sys
import os
import time
import logging
import pandas as pd

# visualization
import panel as pn
import bokeh.server.views.ws
import dashboard_parts

# style
file = os.path.join(
    os.path.dirname(__file__),
    "style",
    "dashboard_style.css"
)
with open(file) as f:
    css = f.read()

pn.extension(raw_css=[css])

# global VARIABLES
SERVER = None
TAB_COUNTER = 0


header = dashboard_parts.HeaderWidget(
    title='AlphaQuant',
    img_folder_path=os.path.join(
        os.path.dirname(__file__),
        "img",
    ),
    github_url='https://github.com/MannLabs/alphaquant'
)

def run():
    global SERVER

    LAYOUT = pn.Column(
        header.create(),
        sizing_mode='stretch_width'
    )

    original_open = bokeh.server.views.ws.WSHandler.open
    bokeh.server.views.ws.WSHandler.open = open_browser_tab(original_open)
    original_on_close = bokeh.server.views.ws.WSHandler.on_close
    bokeh.server.views.ws.WSHandler.on_close = close_browser_tab(
        original_on_close
    )

    SERVER = LAYOUT.show(threaded=True, title='AlphaQuant')
    SERVER.join()

def open_browser_tab(func):
    def wrapper(*args, **kwargs):
        global TAB_COUNTER
        TAB_COUNTER += 1
        return func(*args, **kwargs)
    return wrapper


def close_browser_tab(func):
    def wrapper(*args, **kwargs):
        global TAB_COUNTER
        TAB_COUNTER -= 1
        return_value = func(*args, **kwargs)
        if TAB_COUNTER == 0:
            quit_server()
        return return_value
    return wrapper


def quit_server():
    logging.info("Quitting server...")
    SERVER.stop()



if __name__ == '__main__':
    run()
