import matplotlib.pyplot as plt
import pandas as pd
import argparse
import seaborn as sns
import os
from typing import Any, Dict, List, Tuple
from lib.util import *


class CustomPlotSession:
    FIGURE_PATH = os.path.join("figures")

    def __init__(
        self,
        data: Any,
        setname: str,
        x_axis: str,
        y_axis: str,
        subfigures: List[str],
        styles: List[str],
        filters: List[Dict[str, Any]],
        anti_filters: List[Dict[str, Any]],
        additional_lines_columns: List[str],
    ) -> None:
        self.data = data
        self.setname = setname
        self.x_axis = x_axis
        self.y_axis = y_axis
        self.subfigures = subfigures
        self.styles = styles
        self.filters = filters
        self.anti_filters = anti_filters
        self.additional_lines_columns = additional_lines_columns
        self.filters.append({"set_name": setname})
        self.basename = append_date(setname)
        self.dirpath = check_make_dir(os.path.join(self.FIGURE_PATH, self.basename), 0)
        self.limits = []

    def set_limits(self) -> None:
        max_x = self.data[self.x_axis].max()
        min_x = int(self.data[self.x_axis].min())

        max_y = self.data[self.y_axis].max()
        self.limits = (0, 10, 0.8, 1)

    def apply_filters(self) -> None:
        mask = pd.Series([True] * len(self.data))
        for filter_dict in self.filters:
            for key, val in filter_dict.items():
                mask &= self.data[key] == val
        for filter_dict in self.anti_filters:
            for key, val in filter_dict.items():
                mask &= self.data[key] != val
        self.data = self.data[mask]

    def set_seaborn_settings(self):
        sns.set_context("paper")
        # palette = sns.color_palette("deep", as_cmap=True)
        sns.set(style="darkgrid")

    def save_single_figure(
        self, group_data: Any, style: str, index: int, subfigure: List[Any]
    ) -> None:
        plt.ylim(self.limits[2], self.limits[3])
        plt.figure(figsize=(5, 5))
        relplot = sns.relplot(
            x=self.x_axis,
            y=self.y_axis,
            kind="line",
            style=style,
            data=group_data,
        )
        for col in self.additional_lines_columns:
            sns.lineplot(
                x=self.x_axis,
                y=col,  # specify the additional column here
                color="red",  # set the color of the line to red
                ci=None,  # remove the standard deviation error shadow
                data=group_data,
            )
        title = COLLOG["pretty"][self.y_axis] + " vs " + COLLOG["pretty"][self.x_axis]
        plt.title(title)
        plt.xlabel(COLLOG["pretty"][self.x_axis])
        plt.ylabel(COLLOG["pretty"][self.y_axis])
        x_list = list(range(self.limits[0], self.limits[1] + 1))
        plt.xticks(x_list)
        figure_name = (
            self.y_axis
            + "_vs_"
            + self.x_axis
            + "_style-"
            + style
            + "_sub-"
            + "-".join(self.subfigures)
            + "_"
            + "-".join([str(subfigure) for subfigure in subfigure])
            + "_"
            + str(index)
            + ".png"
        )
        figure_path_png = os.path.join(self.dirpath, figure_name)
        plt.savefig(figure_path_png)

    def save_plots(self) -> None:
        self.apply_filters()
        self.set_seaborn_settings()
        self.set_limits()
        print(self.data)
        # TODO: awkward when styles is a singleton
        for style in self.styles:
            grouped = self.data.groupby(self.subfigures)
            index = 0
            for key, group in grouped:
                if not isinstance(key, tuple):
                    key = (key,)
                index += 1
                self.save_single_figure(pd.DataFrame(group), style, index, list(key))


def main():
    # UI/Parser takes arguments from user and holds on to them.
    parser = argparse.ArgumentParser(description="Filter CSV file based on criteria.")
    parser.add_argument("file_path", help="Path to the CSV file")
    parser.add_argument("--set_name", type=str, help="Filter by set name")
    args = parser.parse_args()
    # Parser hands existing file path to Reader, reads and stores existing dataframe
    try:
        df = pd.read_csv(args.file_path)
    except FileNotFoundError:
        print("Error: File not found.")
        return

    # Reader hands raw data to Cleaner (??)
    cleanup_to_processed(df)

    # Cleaner hands clean data to Finisher (??)
    data_df = cleanup_to_finished(df)
    pretty_df = cleanup_to_finished(data_df)
    del df
    del data_df


    # Finisher hands finished data to CustomPlotSession
    session1 = CustomPlotSession(
        pretty_df,
        args.set_name,
        "budget",
        "empirical_suboptimal_ratio",
        ["policies"],
        ["solver"],
        [
            {"solver": "GREEDY"},
            {"subsolver": "MIP"},
            {"k_zero": 3},
            {"scenarios": 5},
            {"nodes": 79},
        ],
        [
            {"policies": 0},
            {"policies": 1},
            {"policies": 4},
            {"policies": 5},
        ],
        ["exact_alpha"],
    )
    session1.save_plots()
    # Running time vs budget - GREEDY, a curve per n, a plot per k
    # session1 = CustomPlotSession(
    #     pretty_df,
    #     args.set_name,
    #     "budget",
    #     "time_s",
    #     ["policies"],
    #     ["nodes"],
    #     [
    #         {"solver": "GREEDY"},
    #         {"k_zero": 5},
    #     ],
    #     [
    #         {"policies": 0},
    #         {"policies": 1},
    #         {"policies": 4},
    #         {"policies": 5},
    #         {"nodes": 202},
    #         {"nodes": 502},
    #     ],
    # )
    # # Running time vs budget - GREEDY, a curve per k, a plot per n
    # session2 = CustomPlotSession(
    #     pretty_df,
    #     args.set_name,
    #     "budget",
    #     "time_s",
    #     ["nodes"],
    #     ["policies"],
    #     [
    #         {"solver": "GREEDY"},
    #         {"k_zero": 5},
    #     ],
    #     [
    #         {"policies": 0},
    #         {"policies": 1},
    #         {"policies": 4},
    #         {"policies": 5},
    #         {"nodes": 202},
    #         {"nodes": 502},
    #     ],
    # )
    # session1.save_plots()
    # session2.save_plots()
    # Approximation ratio vs budget - GREEDY algorithm, a curve per subsolver, a plot per k
    # session2 = CustomPlotSession()

    # Running time vs budget


if __name__ == "__main__":
    main()
