#!/usr/bin/env python3

# Main file for generating polynomials from the depth map

# $ python im.py [path to IMAGE]

import numpy as np
from polynomial import Polynomial
from PIL import Image
import sys
import os
import tqdm
from pathlib import Path
from multiprocessing.pool import Pool
import argparse


def gethots(n, inds):
    """Create a one hot vector"""
    a = np.zeros(n, dtype=np.uint16)
    a[tuple(i for i in inds)] += 1
    return tuple(a.flatten())


def loop_body(ip):
    """Compute the polynomial for each r"""
    i_j, u, v, arr, l, folder_name, verbose = ip
    i, j = i_j
    nv = 3
    # r = ((i, j), (i, j + 1), (i + 1, j))
    denom = (u[i, j + 1] - u[i, j]) * (v[i + 1, j] - v[i, j]) - (
        u[i + 1, j] - u[i, j]
    ) * (v[i, j + 1] - v[i, j])
    c1 = v[i, j + 1] - v[i, j]
    c3 = -(v[i + 1, j] - v[i, j])
    c2 = -c3 - c1
    pr = Polynomial(
        nvars=nv,
        coeffs={
            (1, 0, 0): c1 / denom,
            (0, 1, 0): c2 / denom,
            (0, 0, 1): c3 / denom,
        },
    )
    c1_b = u[i, j + 1] - u[i, j]
    c3_b = -(u[i + 1, j] - u[i, j])
    c2_b = -c3_b - c1_b
    qr = Polynomial(
        nvars=nv,
        coeffs={
            (1, 0, 0): c1_b / denom,
            (0, 1, 0): c2_b / denom,
            (0, 0, 1): c3_b / denom,
        },
    )
    t = pr * l[0] + qr * l[1] + 1
    g = (pr * pr + qr * qr + 1) * (arr[i, j] ** 2) - t * t
    op = g * g
    op.to_file(f"{folder_name}/{i}_{j}.txt", verbose, i, j)
    return op


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("file_name", type=Path)
    parser.add_argument("-v", action="store_true")
    args = parser.parse_args()
    file_name = args.file_name
    verbose = args.v

    l = [0, 0, 1]

    folder_name = file_name.name.split(".")[0]
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)
    # u and v are coords:
    im = Image.open(file_name)
    arr = np.asarray(im)
    n = arr.shape
    coords_mesh = np.dstack(
        np.meshgrid(
            np.arange(n[0], dtype=np.float32), np.arange(n[1], dtype=np.float32)
        )
    )
    u = coords_mesh[:, :, 0]
    v = coords_mesh[:, :, 1]
    with Pool() as p:
        for i in tqdm.tqdm(
            p.imap_unordered(
                loop_body,
                [
                    ((i, j), u, v, arr, l, folder_name, verbose)
                    for i in range(n[0] - 1)
                    for j in range(n[1] - 1)
                ],
            ),
            total=(n[0] - 1) * (n[1] - 1),
        ):
            pass


if __name__ == "__main__":
    main()
