#!/usr/bin/env python3

import numpy as np
import math
from PIL import Image

allowed_x = [0,1,2,3,4,6,8,10]
allowed_y = {1.2 : (1040,16), 1 : (16,16), 2 : (16,1040)}
max_width = max(allowed_x)

buffer = np.full((2048,2048,3),128,dtype=np.uint8)
buffer[:,:,2] = 255

unit = 100

def fill(width, height, gutter = 3):
	output = np.zeros((int(height*unit),int(width*unit),3))
	if width and height:
		xBuf = np.zeros(int(width*unit))
		xBuf[:gutter] = np.arange(-gutter,0)
		xBuf[-gutter:] = np.arange(1,1+gutter)
		xBuf /= float(gutter)
		yBuf = np.zeros(int(height*unit))
		yBuf[:gutter] = np.arange(-gutter,0)
		yBuf[-gutter:] = np.arange(1,1+gutter)
		yBuf /= float(gutter)
		
		xv, yv = np.meshgrid(xBuf, yBuf)
		
		coords = np.stack([xv, yv],-1)
		
		#physics convention
		phi = np.where(xv**2 + yv**2 > 0, np.arctan2(yv, xv), 0)
		theta = np.maximum(np.abs(xv),np.abs(yv)) * math.pi / 4
		
		
		output[:,:,0] = np.sin(theta)*np.cos(phi)
		output[:,:,1] = np.sin(theta)*np.sin(phi)
		output[:,:,2] = np.cos(theta)
		
	return (np.round(output * 127) + 128).astype(np.uint8)
	
for (height, (startRow, startCol)) in allowed_y.items():
	for brickRow in range(int(len(allowed_x) / 2)):
		smallWidth = allowed_x[brickRow]
		bigWidth = allowed_x[-(1+brickRow)]
		
		rowSpan = slice(int(startRow + height * unit * brickRow), int(startRow + height * unit * (1 + brickRow)))
		smallColSpan = slice(startCol,startCol + smallWidth * unit)
		bigColSpan = slice(startCol + (max_width - bigWidth) * unit,startCol + max_width * unit)
		buffer[rowSpan,smallColSpan,:] = fill(smallWidth,height)
		buffer[rowSpan,bigColSpan,:] = fill(bigWidth,height)

buffer = np.flip(buffer,0)
brickNorms = Image.fromarray(buffer)
brickNorms.save("/Users/thomas/Downloads/brickNorms.png")