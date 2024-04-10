import matplotlib.pyplot as plt
import pandas as pd
import argparse
import seaborn as sns
import os
from typing import Any, Dict, List, Tuple
from lib.util import (
    append_date,
    check_make_dir,
    common_cleanup,
    get_solver_from_flag,
    get_subsolver_from_flag,
    COLLOG,
)


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
    ) -> None:
        self.data = data
        self.setname = setname
        self.x_axis = x_axis
        self.y_axis = y_axis
        self.subfigures = subfigures
        self.styles = styles
        self.filters = filters
        self.anti_filters = anti_filters
        self.filters.append({"set_name": setname})
        self.basename = append_date(setname)
        self.dirpath = check_make_dir(os.path.join(self.FIGURE_PATH, self.basename), 0)

    def set_limits(self) -> None:
        max_x = self.data[self.x_axis].max()
        min_x = int(self.data[self.x_axis].min())

        max_y = self.data[self.y_axis].max()
        # self.limits = (min_x, max_x, 0, max_y)
        self.limits = (0, 10, 0, 50)

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
        x_list = list(range(self.limits[0], self.limits[1] + 1))
        plt.figure(figsize=(5, 5))
        relplot = sns.relplot(
            x=self.x_axis,
            y=self.y_axis,
            kind="line",
            style=style,
            data=group_data,
        )

        title = COLLOG["pretty"][self.y_axis] + " vs " + COLLOG["pretty"][self.x_axis]
        plt.title(title)
        plt.xlabel(COLLOG["pretty"][self.x_axis])
        plt.ylabel(COLLOG["pretty"][self.y_axis])
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
        for style in self.styles:
            grouped = self.data.groupby(self.subfigures)
            index = 0
            for key, group in grouped:
                if not isinstance(key, tuple):
                    key = (key,)
                index += 1
                self.save_single_figure(pd.DataFrame(group), style, index, list(key))


def main():
    parser = argparse.ArgumentParser(description="Filter CSV file based on criteria.")
    parser.add_argument("file_path", help="Path to the CSV file")
    parser.add_argument("--set_name", type=str, help="Filter by set name")
    args = parser.parse_args()
    try:
        df = pd.read_csv(args.file_path)
    except FileNotFoundError:
        print("Error: File not found.")
        return

    # COMMON CLEANUP
    common_cleanup(df)
    # Running time vs budget - GREEDY, a curve per n, a plot per k
    session1 = CustomPlotSession(
        df,
        args.set_name,
        "budget",
        "time_s",
        ["policies"],
        ["nodes"],
        [
            {"solver": "GREEDY"},
            {"k_zero": 5},
        ],
        [
            {"policies": 0},
            {"policies": 1},
            {"policies": 4},
            {"policies": 5},
        ],
    )
    # Running time vs budget - GREEDY, a curve per k, a plot per n
    session2 = CustomPlotSession(
        df,
        args.set_name,
        "budget",
        "time_s",
        ["nodes"],
        ["policies"],
        [
            {"solver": "GREEDY"},
            {"k_zero": 5},
        ],
        [
            {"policies": 0},
            {"policies": 1},
            {"policies": 4},
            {"policies": 5},
        ],
    )
    session1.save_plots()
    session2.save_plots()
    # Approximation ratio vs budget - GREEDY algorithm, a curve per subsolver, a plot per k
    # session2 = CustomPlotSession()

    # Running time vs budget


if __name__ == "__main__":
    main()
