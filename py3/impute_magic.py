#!/usr/bin/env python
import os

import magic
import pandas as pd
import matplotlib.pyplot as plt

projdir = ""
test_data = os.path.join(projdir, "temps", "test_data.csv")
x_raw = pd.read_csv(test_data)

magic_ops = magic.MAGIC()
x_magic = magic_ops.fig_transform(x_raw, genes=["VIM"])
