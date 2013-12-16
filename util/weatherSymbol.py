from math import atan, cos, sin, sqrt
import re

import numpy as np

from matplotlib.lines import Line2D
from matplotlib.patches import Circle
from matplotlib.patches import Wedge
from matplotlib.patches import Rectangle
from matplotlib.patches import Polygon
from matplotlib.patches import Arc
import pylab

from mpl_toolkits.basemap import Basemap

line_width = 1
pi = 4 * atan(1)
deg2rad = pi / 180

def _vec_add(A, B):
    result = []
    for idx in range(len(A)):
        result.append(A[idx] + B[idx])
    return result

def drawRainSymbol(center, rad, rotation=0.0, color='k'):
    pylab.gca().add_patch(Circle(center, radius=rad, edgecolor=color, facecolor=color, linewidth=.5 * line_width))
    return

def drawDrizzleSymbol(center, rad, rotation=0.0, color='k'):
    # Proportion parameter for circle vs. triangle
    alpha = 1. / 3. 

    # x and y coordinates of the point of tangency between one leg of the triangle and the circle.
    tangent_x = rad * (2 * sqrt(alpha) * (1 - alpha)) / (1 + alpha)
    tangent_y = rad * (alpha - (1 - alpha) ** 2 / (1 + alpha))

    # Point list for the triangle
    poly_points = np.array([[0, -rad], [0, rad * alpha], [tangent_x, tangent_y]]) + np.array(center).reshape((1, 2)).repeat(3, axis=0)

    pylab.gca().add_patch(Circle(_vec_add(center, (0, rad * alpha)), radius=rad * (1 - alpha), edgecolor=color, facecolor=color, linewidth=.5 * line_width))
    pylab.gca().add_patch(Polygon(poly_points, edgecolor=color, facecolor=color, linewidth=.5 * line_width))
    return

def drawSnowSymbol(center, rad, rotation=0.0, color='k'):
    x, y = center
    for theta in [ 0, 60, 120 ]:
        start_x = x - rad * cos((rotation + theta) * deg2rad)
        end_x = x + rad * cos((rotation + theta) * deg2rad)

        start_y = y - rad * sin((rotation + theta) * deg2rad)
        end_y = y + rad * sin((rotation + theta) * deg2rad)

        pylab.gca().add_line(Line2D([start_x, end_x], [start_y, end_y], color=color, linewidth=.5 * line_width, solid_capstyle='butt'))
    return

def drawThunderSymbol(center, rad, heavy=False, rotation=0.0, color='k'):
    x, y = center
    x_pts = []
    y_pts = []
    for theta in [ 225, 135, 45 ]:
        x_pts.append(x + rad * cos((rotation + theta) * deg2rad))
        y_pts.append(y + rad * sin((rotation + theta) * deg2rad))

    if False: #heavy:
#       x_pts.extend([x + rad * cos((rotation + 90) * deg2rad) / 3, x + rad * cos((rotation + 90) * degrad), x - ])
#       y_pts.extend([y + rad / 3., y - rad / 3., y-rad])
        pylab.gca().add_line(Line2D(x_pts, y_pts, color=color, linewidth=.5*line_width, solid_capstyle='butt'))
    else:
        x_pts.extend([x, x + rad * cos((rotation - 45) * deg2rad)])
        y_pts.extend([y, y + rad * sin((rotation - 45) * deg2rad)])
        pylab.gca().add_line(Line2D(x_pts, y_pts, color=color, linewidth=.5*line_width, solid_capstyle='butt'))
    return

def drawPelletSymbol(center, rad, sleet=False, rotation=0.0, color='k'):
    x, y = center

    angles = np.array([ 90, 210, 330 ])
    poly_points = np.zeros((len(angles), 2))

    poly_points[:,0] = x + rad * np.cos((angles + rotation) * deg2rad)
    poly_points[:,1] = y + rad * np.sin((angles + rotation) * deg2rad)

    pylab.gca().add_patch(Polygon(poly_points, edgecolor=color, facecolor='w', linewidth=.5 * line_width))
    if sleet:
        drawRainSymbol(center, .25 * rad)
    return

def drawRSDPrecip(center, intensity, rad, rain=False, drizzle=False, snow=False, intermittent=False, color='k'):
    symbol_fn = { 'rain':drawRainSymbol, 'drizzle':drawDrizzleSymbol, 'snow':drawSnowSymbol }

    symbol_rad = rad * (1 - sin(45 * deg2rad))

    if (rain ^ drizzle ^ snow) and not (rain and drizzle and snow):
        if   rain:    symbol = 'rain'
        elif snow:    symbol = 'snow'
        elif drizzle: symbol = 'drizzle'

        if intensity == '-':
            if intermittent:
                symbol_fn[symbol](center, rad, color=color)
            else:
                symbol_fn[symbol](_vec_add(center, (-0.75 * rad, 0)), symbol_rad, color=color)
                symbol_fn[symbol](_vec_add(center, ( 0.75 * rad, 0)), symbol_rad, color=color)
        elif intensity == '':
            if intermittent:
                symbol_fn[symbol](_vec_add(center, (0, -0.5 * rad)), symbol_rad, color=color)
                symbol_fn[symbol](_vec_add(center, (0,  0.5 * rad)), symbol_rad, color=color)
            else:
                symbol_fn[symbol](_vec_add(center, (-0.75 * rad, -0.5 * rad)), symbol_rad, color=color)
                symbol_fn[symbol](_vec_add(center, ( 0.75 * rad, -0.5 *rad)), symbol_rad, color=color)
                symbol_fn[symbol](_vec_add(center, (          0,  0.5 * rad)), symbol_rad, color=color)
        elif intensity == '+':
            if intermittent:
                symbol_fn[symbol](_vec_add(center, (0, -rad)), symbol_rad, color=color)
                symbol_fn[symbol](center, symbol_rad, color=color)
                symbol_fn[symbol](_vec_add(center, (0,  rad)), symbol_rad, color=color)
            else:
                symbol_fn[symbol](_vec_add(center, (          0, -rad)), symbol_rad, color=color)
                symbol_fn[symbol](_vec_add(center, (-0.75 * rad,  0  )), symbol_rad, color=color)
                symbol_fn[symbol](_vec_add(center, ( 0.75 * rad,  0  )), symbol_rad, color=color)
                symbol_fn[symbol](_vec_add(center, (          0,  rad)), symbol_rad, color=color)
    elif (rain and drizzle) or (rain and snow):
        if   rain and snow:    symbol = 'snow'
        elif rain and drizzle: symbol = 'drizzle'

        if intensity == '-':
                symbol_fn[symbol](_vec_add(center, (0, -0.5 * rad)), symbol_rad, color=color)
                drawRainSymbol   (_vec_add(center, (0,  0.5 * rad)), symbol_rad, color=color)
        else:
                symbol_fn[symbol](_vec_add(center, (0, -rad)), symbol_rad, color=color)
                drawRainSymbol   (center, rad)
                symbol_fn[symbol](_vec_add(center, (0,  rad)), symbol_rad, color=color)
    else:
        symbol_list = []
        if rain:    symbol_list.append('rain')
        if snow:    symbol_list.append('snow')
        if drizzle: symbol_list.append('drizzle')
        print "Warning: Symbol combination (%s) not recognized.  Nothing drawn." % ", ".join(symbol_list)
    return

def drawFreezingPrecip(center, intensity, symbol, rad, color='k'):
    symbol_fn = { 'rain':drawRainSymbol, 'drizzle':drawDrizzleSymbol }
    pylab.gca().add_patch(Arc(_vec_add(center, (-.25 * rad, 0)), .5 * rad, .5 * rad, 0., 0., 180., color=color, linewidth=.5 * line_width))
    pylab.gca().add_patch(Arc(_vec_add(center, ( .25 * rad, 0)), .5 * rad, .5 * rad, 0., 180., 360., color=color, linewidth=.5 * line_width))

    symbol_fn[symbol](_vec_add(center, (-.25 * rad, 0)), .09 * rad, color=color)
    if intensity != '-':
        symbol_fn[symbol](_vec_add(center, ( .25 * rad, 0)), .09 * rad, color=color)
    return

def drawThunderPrecip(center, intensity, symbol, rad, color='k'):
    symbol_fn = { 'rain':drawRainSymbol, 'snow':drawSnowSymbol, 'hail':lambda c, r: drawPelletSymbol(c, r, False), 'sleet':lambda c, r: drawPelletSymbol(c, r, True) }
    drawThunderSymbol(center, rad, heavy=intensity == '+', color=color)
    symbol_fn[symbol](_vec_add(center, (0, rad * (1 + sin(45 * deg2rad) / 2))), rad * (1 - sin(45 * deg2rad)), color=color) 
    return

def drawSkyObscuration(center, rad, mist=False, freezing=False, color='k'):
    x, y = center

    if mist: num_lines = 2
    else:    num_lines = 3

    for idx in range(num_lines):
        x_points = [x - rad,                                x + rad                               ]
        y_points = [y + rad * (idx - .5 * (num_lines - 1)), y + rad * (idx - .5 * (num_lines - 1))]
        pylab.gca().add_line(Line2D(x_points, y_points, linewidth=.5 * line_width, color=color, solid_capstyle='butt'))

    if freezing and not mist:
        triangle_width = .75
        x_points = [x - triangle_width * rad, x,       x + triangle_width * rad]
        y_points = [y + rad,                  y - rad, y + rad                 ]
        pylab.gca().add_line(Line2D(x_points, y_points, linewidth=.5 * line_width, color=color, solid_capstyle='butt'))

        poly_points = np.zeros((3,2))
        poly_points[:,0] = [x - triangle_width * rad / 2, x,       x + triangle_width * rad / 2]
        poly_points[:,1] = [y,                            y - rad, y                           ]
        pylab.gca().add_patch(Polygon(poly_points, edgecolor=color, facecolor=color, linewidth=.5 * line_width))
    return

def parseAndPlot(string, center, rad, color='k'):
    if '+' in string:
        qualifier = '+'
    elif '-' in string:
        qualifier = '-'
    else:
        qualifier = ''

    if 'TS' in string:
        if 'RA' in string:
            drawThunderPrecip(center, qualifier, 'rain', rad, color=color)
        elif 'SN' in string:
            drawThunderPrecip(center, qualifier, 'snow', rad, color=color)
        elif 'PL' in string:
            drawThunderPrecip(center, qualifier, 'hail', rad, color=color)
        else:
            drawThunderSymbol(center, rad, heavy=qualifier == '+', color=color)
    elif 'RA' in string:
        if 'FZ' in string:
            drawFreezingPrecip(center, qualifier, 'rain', rad, color=color)
        elif 'SN' in string:
            drawRSDPrecip(center, qualifier, rad, rain=True, snow=True, color=color)
        else:
            drawRSDPrecip(center, qualifier, rad, rain=True, color=color)
    elif 'SN' in string:
        drawRSDPrecip(center, qualifier, rad, snow=True, color=color)
    elif 'PL' in string:
        drawPelletSymbol(center, rad, sleet=True, color=color)
    elif 'DZ' in string:
        if 'FZ' in string:
            drawFreezingPrecip(center, qualifier, 'drizzle', rad, color=color)
        else:
            drawRSDPrecip(center, qualifier, rad, drizzle=True, color=color)
    elif 'FG' in string:
        if 'FZ' in string:
            drawSkyObscuration(center, rad, freezing=True, color=color)
        else:
            drawSkyObscuration(center, rad, color=color)
    elif 'BR' in string:
        if 'FZ' in string:
            drawSkyObscuration(center, rad, mist=True, freezing=True, color=color)
        else:
            drawSkyObscuration(center, rad, mist=True, color=color)
    return

def drawSkyCover(center, rad, oktas=0):
    if oktas < 0 or oktas > 8: return

    x, y = center

    if oktas in [ 7, 8 ]:
        pylab.gca().add_patch(Circle(center, radius=rad, edgecolor='k', facecolor='k', linewidth=line_width))
    else:
        pylab.gca().add_patch(Circle(center, radius=rad, edgecolor='k', facecolor='w', linewidth=line_width))

    if oktas == 7:
        pylab.gca().add_patch(Rectangle((x - .25 * rad, y - .75 * rad), .5 * rad, 1.5 * rad, edgecolor='k', facecolor='w', linewidth=line_width))

    if oktas in [ 2, 3 ]:
        pylab.gca().add_patch(Wedge(center, rad, 0.0, 90.0, edgecolor='k', facecolor='k', linewidth=0.25 * line_width))
    elif oktas in [ 4, 5 ]:
        pylab.gca().add_patch(Wedge(center, rad, -90.0, 90.0, edgecolor='k', facecolor='k', linewidth=0.25 * line_width))
    elif oktas == 6:
        pylab.gca().add_patch(Wedge(center, rad, -180.0, 90.0, edgecolor='k', facecolor='k', linewidth=0.25 * line_width))

    if oktas == 1:
        pylab.gca().add_line(Line2D([x - 0.1 * rad, x - 0.1 * rad], [y + .75 * rad, y - .75 * rad], color='k', linewidth=line_width, solid_capstyle='butt'))
    elif oktas == 3:
        pylab.gca().add_line(Line2D([x - 0.1 * rad, x - 0.1 * rad], [y + rad, y - rad], color='k', linewidth=line_width, solid_capstyle='butt'))
    elif oktas == 5:
        pylab.gca().add_line(Line2D([x + rad, x - rad], [y, y], color='k', linewidth=line_width, solid_capstyle='butt'))
    return

if __name__ == "__main__":
    l_scale = .01
    oktas = 2
    pylab.figure(figsize=(12, 8))
    pylab.axes((0, 0, 1, 1))

    map = Basemap(resolution='i', projection='stere', lat_ts=60., 
        lat_0=60., lon_0=-97.5,
        llcrnrlat=30.0, llcrnrlon=-106.0,
        urcrnrlat=40.0, urcrnrlon=-93.0)
    unit_scale = l_scale * max(map.urcrnry - map.llcrnry, map.urcrnrx - map.llcrnrx)

    map.drawcoastlines()
    map.drawcountries()
    map.drawstates()

    point = map(-97.4167, 35.2167)

    drawSkyCover(point, .7 * unit_scale, oktas)
#   drawRSDPrecip(_vec_add(point, (-1.5 * unit_scale, 0)), '+', .3 * unit_scale, snow=True, rain=True, intermittent=False)
#   drawFreezingPrecip(_vec_add(point, (-2 * unit_scale, 0)), '-', 'drizzle', 1.5 * unit_scale)
#   drawThunderPrecip(_vec_add(point, (-1.5 * unit_scale, 0)), '-', 'hail', .8 * unit_scale)
#   drawPelletSymbol(_vec_add(point, (-1.5 * unit_scale, 0)), .8 * unit_scale, sleet=True)
    drawSkyObscuration(_vec_add(point, (-1.5 * unit_scale, 0)), .5 * unit_scale, freezing=True)
    pylab.barbs(point[0], point[1], 5., 10.)

    pylab.savefig("mpl_test.png")
