import sys
import os
import time
import logging
import pandas as pd
import param

# visualization
import panel as pn
import bokeh.server.views.ws
import alphaquant.ui.dashboard_parts_run_pipeline as dashboard_parts
import alphaquant.ui.gui_textfields as gui_textfields
import alphaquant.ui.dashboad_parts_plots_basic as dashboad_parts_plots_basic
import alphaquant.ui.dashboard_parts_plots_proteoforms as dashboard_parts_plots_proteoforms


def get_css_style(
    file_name="dashboard_style.css",
    directory=os.path.join(
        os.path.dirname(__file__),
        "..",
        "resources",
        "style"
    )
):
    file = os.path.join(
        directory,
        file_name
    )
    with open(file) as f:
        return f.read()


def init_panel():
    pn.extension(raw_css=[get_css_style()])
    pn.extension('plotly')


# style
init_panel()


class DashboardState(param.Parameterized):
    """Central state manager for the dashboard."""
    results_dir = param.String(default="")
    samplemap_df = param.DataFrame(default=pd.DataFrame())
    analysis_file = param.String(default="")
    condition_pairs = param.List(default=[])
    selected_condition_pair = param.String(default="")

    def __init__(self, **params):
        super().__init__(**params)
        self.subscribers = []

    def add_subscriber(self, subscriber):
        """Add a tab/component that should be notified of state changes."""
        self.subscribers.append(subscriber)

    def notify_subscribers(self, changed_param, value=None):
        """Notify all subscribers of a state change."""
        print(f"\n=== State Manager: Notifying Subscribers ===")
        print(f"Changed parameter: {changed_param}")
        print(f"Value type: {type(value)}")
        print(f"Value: {value}")
        if not hasattr(self, '_subscribers'):
            print("No subscribers registered")
            return

        print(f"Number of subscribers: {len(self._subscribers)}")
        for subscriber in self._subscribers:
            print(f"Notifying subscriber: {type(subscriber).__name__}")
            try:
                getattr(subscriber, f'on_{changed_param}_changed')(value)
                print(f"Successfully notified {type(subscriber).__name__}")
            except Exception as e:
                print(f"Error notifying {type(subscriber).__name__}: {str(e)}")
        print("=== Finished Notifying Subscribers ===\n")


class GUI(object):
    # TODO: import from alphabase

    def __init__(
        self,
        name,
        github_url,
        run_in_background=False,
        automatic_close=True,
    ):
        self.name = name
        self.tab_counter = 0
        self.header = dashboard_parts.HeaderWidget(
            name,
            os.path.join(
                os.path.dirname(__file__),
                "..","resources",
                "img",
            ),
            github_url
        )
        self.layout = pn.Column(
            self.header.create(),
            sizing_mode='stretch_width',
            min_width=1270
        )
        self.run_in_background = run_in_background
        self.automatic_close = automatic_close

    def start_server(self, run_in_background=False):
        if self.automatic_close:
            bokeh_ws_handler = bokeh.server.views.ws.WSHandler
            self.bokeh_server_open = bokeh_ws_handler.open
            bokeh_ws_handler.open = self.__open_browser_tab(
                self.bokeh_server_open
            )
            self.bokeh_server_on_close = bokeh_ws_handler.on_close
            bokeh_ws_handler.on_close = self.__close_browser_tab(
                self.bokeh_server_on_close
            )
        self.server = self.layout.show(threaded=True, title=self.name)
        if not run_in_background:
            self.server.join()
        elif not self.run_in_background:
            self.server.join()

    def __open_browser_tab(self, func):
        def wrapper(*args, **kwargs):
            self.tab_counter += 1
            return func(*args, **kwargs)
        return wrapper

    def __close_browser_tab(self, func):
        def wrapper(*args, **kwargs):
            self.tab_counter -= 1
            return_value = func(*args, **kwargs)
            if self.tab_counter == 0:
                self.stop_server()
            return return_value
        return wrapper

    def stop_server(self):
        logging.info("Stopping server...")
        self.server.stop()
        if self.automatic_close:
            bokeh_ws_handler = bokeh.server.views.ws.WSHandler
            bokeh_ws_handler.open = self.bokeh_server_open
            bokeh_ws_handler.on_close = self.bokeh_server_on_close


class AlphaQuantGUI(GUI):
    def __init__(self, start_server=False):
        super().__init__(
            name="AlphaQuant",
            github_url='https://github.com/MannLabs/alphaquant',
        )

        # Create central state manager
        self.state = DashboardState()

        # Create components/tabs
        self.run_pipeline = dashboard_parts.RunPipeline(state=self.state)
        self.plotting_tab = dashboad_parts_plots_basic.PlottingTab(state=self.state)
        self.proteoform_tab = dashboard_parts_plots_proteoforms.ProteoformPlottingTab(state=self.state)

        # Register components as subscribers
        self.state.add_subscriber(self.run_pipeline)
        self.state.add_subscriber(self.plotting_tab)
        self.state.add_subscriber(self.proteoform_tab)

        self.project_description = """<div style="color: #2F4F4F; font-size: 1.3em; margin-top: -10px; margin-bottom: 20px;">AlphaQuant is an open-source package for peptide-resolved detection of protein abundance changes.</div>"""

        # Create a centered row for the project description
        self.description_row = pn.Row(
            pn.Spacer(sizing_mode='stretch_width'),
            pn.pane.HTML(self.project_description, align='center'),
            pn.Spacer(sizing_mode='stretch_width'),
            sizing_mode='stretch_width',
            margin=(0, 0, 20, 0)  # top, right, bottom, left
        )

        # Create instructions card
        self.instructions_card = pn.Card(
            "### Instructions",
            gui_textfields.Descriptions.project_instruction,
            gui_textfields.Cards.spectronaut,
            gui_textfields.Cards.diann,
            gui_textfields.Cards.alphapept,
            gui_textfields.Cards.maxquant,
            title='Instructions',
            collapsed=True,
            margin=(5, 5, 5, 5),
            sizing_mode='fixed',
        )

        # Wrap instructions card in a Row for horizontal centering
        self.instructions_row = pn.Row(
            pn.Spacer(sizing_mode='stretch_width'),
            self.instructions_card,
            pn.Spacer(sizing_mode='stretch_width'),
            sizing_mode='stretch_width'
        )

        # ERROR/WARNING MESSAGES
        self.error_message_upload = "The selected file can't be uploaded. Please check the instructions for data uploading."

        # Create initial tabs
        self.tab_layout = pn.Tabs(
            ('Run Pipeline', self.run_pipeline.create()),
            ('Basic Plots', self.plotting_tab.panel()),
            ('Proteoform Plots', self.proteoform_tab.panel()),
            dynamic=True,
            tabs_location='above',
            sizing_mode='stretch_width'
        )

        self.layout += [
            self.description_row,
            self.instructions_row,
            self.tab_layout
        ]

        # Add process management and cleanup
        self._server = None
        self._workers = []

        if start_server:
            self.start_server()

    def start_server(self):
        try:
            # Configure process pool with explicit resource limits
            import multiprocessing as mp
            self._workers = mp.Pool(
                processes=min(mp.cpu_count(), 4),  # Limit max processes
                maxtasksperchild=100  # Restart workers periodically
            )

            # Start server with proper cleanup handling
            self._server = pn.serve(
                self.layout,
                port=0,  # Dynamic port allocation
                show=True,
                websocket_origin="*",
                threaded=True
            )

            # Register cleanup handlers
            import atexit
            atexit.register(self.cleanup)

        except Exception as e:
            self.cleanup()
            logging.error(f"Server startup failed: {str(e)}")
            raise

    def cleanup(self):
        # Graceful cleanup of resources
        if self._workers:
            self._workers.terminate()
            self._workers.join()
            self._workers = None

        if self._server:
            self._server.stop()
            self._server = None


def run():
    try:
        # Add error handling wrapper
        gui = AlphaQuantGUI(start_server=True)
        return gui
    except Exception as e:
        logging.error(f"Failed to start GUI: {str(e)}")
        raise


if __name__ == '__main__':
    run()
