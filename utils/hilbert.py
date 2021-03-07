#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
    This file is for generating hilbert space filling curve for arbitrary shape matrix and its inverse mapping.
    Usage:
        h = hilbert.hilbert_get(x, y)
        x: number of columns
        y: number of rows
"""

l = []


def hilbert_get(x, y):
    global l
    x_old = 0
    y_old = 0
    spacefill(1, 1, x, y)
    mesh = l
    l = []
    return mesh


def render(x0, y0, direction):
    global l
    x_new = x0
    y_new = y0
    l.append([int(y_new - 1), int(x_new - 1)])
    x_old = x_new
    y_old = y_new


def spacefill(ll, tt, ww, hh):
    if hh > ww:
        if hh % 2 == 1 and ww % 2 == 0:
            go(ll, tt, ww, 0, 0, hh, "m")
        else:
            go(ll, tt, ww, 0, 0, hh, "r")
    else:
        if ww % 2 == 1 and hh % 2 == 0:
            go(ll, tt, ww, 0, 0, hh, "m")
        else:
            go(ll, tt, ww, 0, 0, hh, "l")


def go(x0, y0, dxl, dyl, dxr, dyr, direction):
    if abs((dxl+dyl)*(dxr+dyr)) <= 6:
        if abs(dxl+dyl) == 1:
            ddx = int(dxr/abs(dxr+dyr))
            ddy = int(dyr/abs(dxr+dyr))
            for ii in range(abs(dxr+dyr)):
                render(x0+ii*ddx+(dxl+ddx-1)/2, y0+ii*ddy+(dyl+ddy-1)/2, direction)
            return False
        if abs(dxr+dyr) == 1:
            ddx = int(dxl/abs(dxl+dyl))
            ddy = int(dyl/abs(dxl+dyl))
            for ii in range(abs(dxl+dyl)):
                render(x0+ii*ddx+(dxr+ddx-1)/2, y0+ii*ddy+(dyr+ddy-1)/2, direction)        
            return False
        if direction == "l":
            ddx = int(dxr/abs(dxr+dyr))
            ddy = int(dyr/abs(dxr+dyr))
            for ii in range(abs(dxr+dyr)):
                render(x0+ii*ddx+(dxl/2+ddx-1)/2, y0+ii*ddy+(dyl/2+ddy-1)/2, direction)
            for ii in range(abs(dxr+dyr)-1, -1, -1):
                render(x0+ii*ddx+(dxl+dxl/2+ddx-1)/2, y0+ii*ddy+(dyl+dyl/2+ddy-1)/2, direction)
            return False
        if direction == "r":
            ddx = int(dxl/abs(dxl+dyl))
            ddy = int(dyl/abs(dxl+dyl))
            for ii in range(abs(dxl+dyl)):
                render(x0+ii*ddx+(dxr/2+ddx-1)/2, y0+ii*ddy+(dyr/2+ddy-1)/2, direction)
            for ii in range(abs(dxl+dyl)-1, -1, -1):
                render(x0+ii*ddx+(dxr+dxr/2+ddx-1)/2, y0+ii*ddy+(dyr+dyr/2+ddy-1)/2, direction)
            return False
        if direction == "m":
            if abs(dxr+dyr) == 3:
                ddx = int(dxr/abs(dxr+dyr))
                ddy = dyr/abs(dxr+dyr)
                render(x0+(dxl/2+ddx-1)/2, y0+(dyl/2+ddy-1)/2, direction)
                render(x0+(dxl+dxl/2+ddx-1)/2, y0+(dyl+dyl/2+ddy-1)/2, direction)
                render(x0+ddx+(dxl+dxl/2+ddx-1)/2, y0+ddy+(dyl+dyl/2+ddy-1)/2, direction)
                render(x0+ddx+(dxl/2+ddx-1)/2, y0+ddy+(dyl/2+ddy-1)/2, direction)
                render(x0+2*ddx+(dxl/2+ddx-1)/2, y0+2*ddy+(dyl/2+ddy-1)/2, direction)
                render(x0+2*ddx+(dxl+dxl/2+ddx-1)/2, y0+2*ddy+(dyl+dyl/2+ddy-1)/2, direction)
                return False
            if abs(dxl+dyl) == 3:
                ddx = dxl/abs(dxl+dyl)
                ddy = dyl/abs(dxl+dyl)
                render(x0+(dxr/2+ddx-1)/2, y0+(dyr/2+ddy-1)/2, direction)
                render(x0+(dxr+dxr/2+ddx-1)/2, y0+(dyr+dyr/2+ddy-1)/2, direction)
                render(x0+ddx+(dxr+dxr/2+ddx-1)/2, y0+ddy+(dyr+dyr/2+ddy-1)/2, direction)
                render(x0+ddx+(dxr/2+ddx-1)/2, y0+ddy+(dyr/2+ddy-1)/2, direction)
                render(x0+2*ddx+(dxr/2+ddx-1)/2, y0+2*ddy+(dyr/2+ddy-1)/2, direction)
                render(x0+2*ddx+(dxr+dxr/2+ddx-1)/2, y0+2*ddy+(dyr+dyr/2+ddy-1)/2, direction)
                return False
        return False
    if 2*(abs(dxl)+abs(dyl)) > 3*(abs(dxr)+abs(dyr)):
        dxl2 = round(dxl/2)
        dyl2 = round(dyl/2)
        if (abs(dxr)+abs(dyr)) % 2 == 0:
            if (abs(dxl)+abs(dyl)) % 2 == 0:
                if direction == "l":
                    if (abs(dxl)+abs(dyl)) % 4 == 0:
                        go(x0, y0, dxl2, dyl2, dxr, dyr, "l")
                        go(x0+dxl2, y0+dyl2, dxl-dxl2, dyl-dyl2, dxr, dyr, "l")
                    else:
                        go(x0, y0, dxl2, dyl2, dxr, dyr, "m")
                        go(x0+dxl2+dxr, y0+dyl2+dyr, -dxr, -dyr, dxl-dxl2, dyl-dyl2, "m") 
                    return False
            else:
                if direction == "m":
                    if (abs(dxl2)+abs(dyl2)) % 2 == 0:
                        go(x0, y0, dxl2, dyl2, dxr, dyr, "l")
                        go(x0+dxl2, y0+dyl2, dxl-dxl2, dyl-dyl2, dxr, dyr, "m")
                    else:
                        go(x0, y0, dxl2, dyl2, dxr, dyr, "m")
                        go(x0+dxl2+dxr, y0+dyl2+dyr, -dxr, -dyr, dxl-dxl2, dyl-dyl2, "r")
                    return False
        else:
            if direction == "l":
                go(x0, y0, dxl2, dyl2, dxr, dyr, "l")
                go(x0+dxl2, y0+dyl2, dxl-dxl2, dyl-dyl2, dxr, dyr, "l")
                return False
            if direction == "m":
                go(x0, y0, dxl2, dyl2, dxr, dyr, "l")
                go(x0+dxl2, y0+dyl2, dxl-dxl2, dyl-dyl2, dxr, dyr, "m")
                return False
    if 2*(abs(dxr)+abs(dyr)) > 3*(abs(dxl)+abs(dyl)):
        dxr2 = round(dxr/2)
        dyr2 = round(dyr/2)
        if (abs(dxl)+abs(dyl)) % 2 == 0:
            if (abs(dxr)+abs(dyr)) % 2 == 0:
                if direction == "r":
                    if (abs(dxr)+abs(dyr)) % 4 == 0:
                        go(x0, y0, dxl, dyl, dxr2, dyr2, "r")
                        go(x0+dxr2, y0+dyr2, dxl, dyl, dxr-dxr2, dyr-dyr2, "r")
                    else:
                        go(x0, y0, dxl, dyl, dxr2, dyr2, "m")
                        go(x0+dxr2+dxl, y0+dyr2+dyl, dxr-dxr2, dyr-dyr2, -dxl, -dyl, "m")
                    return False
            else:
                if direction == "m":
                    if (abs(dxr2)+abs(dyr2)) % 2 == 0:
                        go(x0, y0, dxl, dyl, dxr2, dyr2, "r")
                        go(x0+dxr2, y0+dyr2, dxl, dyl, dxr-dxr2, dyr-dyr2, "m")
                    else:
                        go(x0, y0, dxl, dyl, dxr2, dyr2, "m")
                        go(x0+dxr2+dxl, y0+dyr2+dyl, dxr-dxr2, dyr-dyr2, -dxl, -dyl, "l")
                    return False
        else:
            if direction == "r":
                go(x0, y0, dxl, dyl, dxr2, dyr2, "r")
                go(x0+dxr2, y0+dyr2, dxl, dyl, dxr-dxr2, dyr-dyr2, "r")
                return False
            if direction == "m":
                go(x0, y0, dxl, dyl, dxr2, dyr2, "r")
                go(x0+dxr2, y0+dyr2, dxl, dyl, dxr-dxr2, dyr-dyr2, "m")
                return False
    if direction == "l" or direction == "r":
        dxl2 = round(dxl/2)
        dyl2 = round(dyl/2)
        dxr2 = round(dxr/2)
        dyr2 = round(dyr/2)
        if abs(dxl+dyl) % 2 == 0 and abs(dxr+dyr) % 2 == 0:
            if abs(dxl2+dyl2+dxr2+dyr2) % 2 == 0:
                if direction == "l":
                    go(x0, y0, dxl2, dyl2, dxr2, dyr2, "r")
                    go(x0+dxr2, y0+dyr2, dxl2, dyl2, dxr-dxr2, dyr-dyr2, "l")
                    go(x0+dxr2+dxl2, y0+dyr2+dyl2, dxl-dxl2, dyl-dyl2, dxr-dxr2, dyr-dyr2, "l")
                    go(x0+dxr2+dxl, y0+dyr2+dyl, dxl2-dxl, dyl2-dyl, -dxr2, -dyr2, "r")
                else:
                    go(x0, y0, dxl2, dyl2, dxr2, dyr2, "l")
                    go(x0+dxl2, y0+dyl2, dxl-dxl2, dyl-dyl2, dxr2, dyr2, "r")
                    go(x0+dxr2+dxl2, y0+dyr2+dyl2, dxl-dxl2, dyl-dyl2, dxr-dxr2, dyr-dyr2, "r")
                    go(x0+dxr+dxl2, y0+dyr+dyl2, -dxl2, -dyl2, dxr2-dxr, dyr2-dyr, "l")
            else:
                if (dxr2+dyr2) % 2 == 0:
                    if direction == "l":
                        go(x0, y0, dxl2, dyl2, dxr2, dyr2, "r")
                        go(x0+dxr2, y0+dyr2, dxl2, dyl2, dxr-dxr2, dyr-dyr2, "m")
                        go(x0+dxr+dxl2, y0+dyr+dyl2, dxr2-dxr, dyr2-dyr, dxl-dxl2, dyl-dyl2, "m")
                        go(x0+dxr2+dxl, y0+dyr2+dyl, dxl2-dxl, dyl2-dyl, -dxr2, -dyr2, "r")
                    else:
                        if dxr2 != 0:
                            dxr2 = dxr2 + 1
                        else:
                            dyr2 = dyr2 + 1 
                        go(x0, y0, dxl2, dyl2, dxr2, dyr2, "l")
                        go(x0+dxl2, y0+dyl2, dxl-dxl2, dyl-dyl2, dxr2, dyr2, "m")
                        go(x0+dxl+dxr2, y0+dyl+dyr2, dxr-dxr2, dyr-dyr2, dxl2-dxl, dyl2-dyl, "m")
                        go(x0+dxl2+dxr, y0+dyl2+dyr, -dxl2, -dyl2, dxr2-dxr, dyr2-dyr, "l")
                else:
                    if direction == "r":
                        go(x0, y0, dxl2, dyl2, dxr2, dyr2, "l")
                        go(x0+dxl2, y0+dyl2, dxl-dxl2, dyl-dyl2, dxr2, dyr2, "m")
                        go(x0+dxl+dxr2, y0+dyl+dyr2, dxr-dxr2, dyr-dyr2, dxl2-dxl, dyl2-dyl, "m")
                        go(x0+dxl2+dxr, y0+dyl2+dyr, -dxl2, -dyl2, dxr2-dxr, dyr2-dyr, "l")
                    else:
                        if dxl2 != 0:
                            dxl2 = dxl2+1
                        else:
                            dyl2 = dyl2+1
                        go(x0, y0, dxl2, dyl2, dxr2, dyr2, "r")
                        go(x0+dxr2, y0+dyr2, dxl2, dyl2, dxr-dxr2, dyr-dyr2, "m")
                        go(x0+dxr+dxl2, y0+dyr+dyl2, dxr2-dxr, dyr2-dyr, dxl-dxl2, dyl-dyl2, "m")
                        go(x0+dxr2+dxl, y0+dyr2+dyl, dxl2-dxl, dyl2-dyl, -dxr2, -dyr2, "r")
        else:
            if abs(dxl+dyl) % 2 != 0 and abs(dxr+dyr) % 2 != 0:
                if dxl2 % 2 != 0:
                    dxl2 = dxl-dxl2
                if dyl2 % 2 != 0:
                    dyl2 = dyl-dyl2
                if dxr2 % 2 != 0:
                    dxr2 = dxr-dxr2
                if dyr2 % 2 != 0:
                    dyr2 = dyr-dyr2
                if direction == "l":
                    go(x0, y0, dxl2, dyl2, dxr2, dyr2, "r")
                    go(x0+dxr2, y0+dyr2, dxl2, dyl2, dxr-dxr2, dyr-dyr2, "m")
                    go(x0+dxr+dxl2, y0+dyr+dyl2, dxr2-dxr, dyr2-dyr, dxl-dxl2, dyl-dyl2, "m")
                    go(x0+dxr2+dxl, y0+dyr2+dyl, dxl2-dxl, dyl2-dyl, -dxr2, -dyr2, "r")
                else:
                    go(x0, y0, dxl2, dyl2, dxr2, dyr2, "l")
                    go(x0+dxl2, y0+dyl2, dxl-dxl2, dyl-dyl2, dxr2, dyr2, "m")
                    go(x0+dxl+dxr2, y0+dyl+dyr2, dxr-dxr2, dyr-dyr2, dxl2-dxl, dyl2-dyl, "m")
                    go(x0+dxl2+dxr, y0+dyl2+dyr, -dxl2, -dyl2, dxr2-dxr, dyr2-dyr, "l")
            else:
                if abs(dxl+dyl) % 2 == 0:
                    if direction == "l":
                        if dxr2 % 2 != 0:
                            dxr2 = dxr-dxr2
                        if dyr2 % 2 != 0:
                            dyr2 = dyr-dyr2
                        if abs(dxl+dyl) > 2:
                            go(x0, y0, dxl2, dyl2, dxr2, dyr2, "r")
                            go(x0+dxr2, y0+dyr2, dxl2, dyl2, dxr-dxr2, dyr-dyr2, "l")
                            go(x0+dxr2+dxl2, y0+dyr2+dyl2, dxl-dxl2, dyl-dyl2, dxr-dxr2, dyr-dyr2, "l")
                            go(x0+dxr2+dxl, y0+dyr2+dyl, dxl2-dxl, dyl2-dyl, -dxr2, -dyr2, "r")
                        else:
                            go(x0, y0, dxl2, dyl2, dxr2, dyr2, "r")
                            go(x0+dxr2, y0+dyr2, dxl2, dyl2, dxr-dxr2, dyr-dyr2, "m")
                            go(x0+dxr+dxl2, y0+dyr+dyl2, dxr2-dxr, dyr2-dyr, dxl-dxl2, dyl-dyl2, "m")
                            go(x0+dxr2+dxl, y0+dyr2+dyl, dxl2-dxl, dyl2-dyl, -dxr2, -dyr2, "r")
                else:
                    if direction == "r":
                        if dxl2 % 2 != 0:
                            dxl2 = dxl-dxl2
                        if dyl2 % 2 != 0:
                            dyl2 = dyl-dyl2
                        if abs(dxr+dyr) > 2:
                            go(x0, y0, dxl2, dyl2, dxr2, dyr2, "l")
                            go(x0+dxl2, y0+dyl2, dxl-dxl2, dyl-dyl2, dxr2, dyr2, "r")
                            go(x0+dxr2+dxl2, y0+dyr2+dyl2, dxl-dxl2, dyl-dyl2, dxr-dxr2, dyr-dyr2, "r")
                            go(x0+dxr+dxl2, y0+dyr+dyl2, -dxl2, -dyl2, dxr2-dxr, dyr2-dyr, "l")
                        else:
                            go(x0, y0, dxl2, dyl2, dxr2, dyr2, "l")
                            go(x0+dxl2, y0+dyl2, dxl-dxl2, dyl-dyl2, dxr2, dyr2, "m")
                            go(x0+dxl+dxr2, y0+dyl+dyr2, dxr-dxr2, dyr-dyr2, dxl2-dxl, dyl2-dyl, "m")
                            go(x0+dxl2+dxr, y0+dyl2+dyr, -dxl2, -dyl2, dxr2-dxr, dyr2-dyr, "l")
    else:
        if abs(dxr+dyr) % 2 == 0:
            dxl2 = round(dxl/3)
            dyl2 = round(dyl/3)
            dxr2 = round(dxr/3)
            dyr2 = round(dyr/3)
            if (dxl2+dyl2) % 2 == 0:
                dxl2 = dxl-2*dxl2
                dyl2 = dyl-2*dyl2
            if (dxr2+dyr2) % 2 == 0:
                if abs(dxr2+dyr2) != 2:
                    if dxr < 0:
                        dxr2 = dxr2+1
                    if dxr > 0:
                        dxr2 = dxr2-1
                    if dyr < 0:
                        dyr2 = dyr2 + 1
                    if dyr > 0:
                        dyr2 = dyr2 - 1
        else:
            dxl2 = round(dxl/3)
            dyl2 = round(dyl/3)
            dxr2 = round(dxr/3)
            dyr2 = round(dyr/3)
            if (dxr2+dyr2) % 2 == 0:
                dxr2 = dxr-2*dxr2
                dyr2 = dyr-2*dyr2
            if (dxl2+dyl2) % 2 == 0:
                if abs(dxl2+dyl2) != 2:
                    if dxl < 0:
                        dxl2 = dxl2 + 1
                    if dxl > 0:
                        dxl2 = dxl2 - 1
                    if dyl < 0:
                        dyl2 = dyl2 + 1
                    if dyl > 0:
                        dyl2 = dyl2 - 1
        if abs(dxl+dyl) < abs(dxr+dyr):
            go(x0, y0, dxl2, dyl2, dxr2, dyr2, "m")
            go(x0+dxl2+dxr2, y0+dyl2+dyr2, -dxr2, -dyr2, dxl-2*dxl2, dyl-2*dyl2, "m")
            go(x0+dxl-dxl2, y0+dyl-dyl2, dxl2, dyl2, dxr2, dyr2, "m")
            go(x0+dxl+dxr2, y0+dyl+dyr2, dxr-2*dxr2, dyr-2*dyr2, -dxl2, -dyl2, "m")
            go(x0+dxr-dxr2+dxl-dxl2, y0+dyr-dyr2+dyl-dyl2, 2*dxl2-dxl, 2*dyl2-dyl, 2*dxr2-dxr, 2*dyr2-dyr, "m")
            go(x0+dxl2+dxr2, y0+dyl2+dyr2, dxr-2*dxr2, dyr-2*dyr2, -dxl2, -dyl2, "m")
            go(x0+dxr-dxr2, y0+dyr-dyr2, dxl2, dyl2, dxr2, dyr2, "m")
            go(x0+dxr+dxl2, y0+dyr+dyl2, -dxr2, -dyr2, dxl-2*dxl2, dyl-2*dyl2, "m")
            go(x0+dxr-dxr2+dxl-dxl2, y0+dyr-dyr2+dyl-dyl2, dxl2, dyl2, dxr2, dyr2, "m")
        else:
            go(x0, y0, dxl2, dyl2, dxr2, dyr2, "m")
            go(x0+dxl2+dxr2, y0+dyl2+dyr2, dxr-2*dxr2, dyr-2*dyr2, -dxl2, -dyl2, "m")
            go(x0+dxr-dxr2, y0+dyr-dyr2, dxl2, dyl2, dxr2, dyr2, "m")
            go(x0+dxr+dxl2, y0+dyr+dyl2, -dxr2, -dyr2, dxl-2*dxl2, dyl-2*dyl2, "m")
            go(x0+dxr-dxr2+dxl-dxl2, y0+dyr-dyr2+dyl-dyl2, 2*dxl2-dxl, 2*dyl2-dyl, 2*dxr2-dxr, 2*dyr2-dyr, "m")
            go(x0+dxl2+dxr2, y0+dyl2+dyr2, -dxr2, -dyr2, dxl-2*dxl2, dyl-2*dyl2, "m")
            go(x0+dxl-dxl2, y0+dyl-dyl2, dxl2, dyl2, dxr2, dyr2, "m")
            go(x0+dxl+dxr2, y0+dyl+dyr2, dxr-2*dxr2, dyr-2*dyr2, -dxl2, -dyl2, "m")
            go(x0+dxr-dxr2+dxl-dxl2, y0+dyr-dyr2+dyl-dyl2, dxl2, dyl2, dxr2, dyr2, "m")
