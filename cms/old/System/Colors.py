import math

'''
Color space manipulations. 

E.g. HSV2RGB, GrayScale, etc... returned in web-safe hex triplets. 
'''

def GrayScale(t):
	return Point2HexColor(0,0,t)

def Point2HexColor(a,lfrac,tfrac):
	[H,S,V] = [math.floor(360*a),lfrac,tfrac]

	RGB = hsvToRGB(H,S,V)

	H = [hex(int(math.floor(255*x))) for x in RGB]

	HEX = [a[a.find('x')+1:] for a in H]
	HEX = ['0' + h if len(h) == 1 else h for h in HEX]
	
	return '#' + ''.join(HEX)

def hsvToRGB(h, s, v):
	"""Convert HSV color space to RGB color space
	
	@param h: Hue
	@param s: Saturation
	@param v: Value
	return (r, g, b)  
	"""
	hi = math.floor(h / 60.0) % 6
	f =  (h / 60.0) - math.floor(h / 60.0)
	p = v * (1.0 - s)
	q = v * (1.0 - (f*s))
	t = v * (1.0 - ((1.0 - f) * s))
	
	D = {0: (v, t, p), 1: (q, v, p), 2: (p, v, t), 3: (p, q, v), 4: (t, p, v), 5: (v, p, q)}
	
	return D[hi]
