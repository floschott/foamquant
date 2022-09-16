def json_rand_dictionary(Ncolors, namecmap, first_color_black=True):
    from matplotlib.colors import LinearSegmentedColormap
    import colorsys
    import numpy as np
    import json
    
    randHSVcolors = [(np.random.uniform(low=0.0, high=1), np.random.uniform(low=0.2, high=1), np.random.uniform(low=0.9, high=1)) for i in range(Ncolors)]
    
    randRGBcolors = []
    for i in range(len(randHSVcolors)):
        RGBcolors = colorsys.hsv_to_rgb(randHSVcolors[i][0], randHSVcolors[i][1], randHSVcolors[i][2])
        
        x0 = i/Ncolors
        if i >0:
            x0 = i/Ncolors#+0.001
        x1 = (i+1)/Ncolors
        
        randRGBcolors.append(x0)
        randRGBcolors.append(RGBcolors[0]); randRGBcolors.append(RGBcolors[1]); randRGBcolors.append(RGBcolors[2])
        randRGBcolors.append(x1)
        randRGBcolors.append(RGBcolors[0]); randRGBcolors.append(RGBcolors[1]); randRGBcolors.append(RGBcolors[2])
    
    #randRGBcolors = np.asarray(randRGBcolors)
    
    
    #print(randRGBcolors)
    if first_color_black == True:
        randRGBcolors[1:4] = [0,0,0]
        randRGBcolors[5:8] = [0,0,0]
    
    
    json_cmap = [
	{
		"ColorSpace" : "HSV",
		"Name" : namecmap,
		"RGBPoints" : randRGBcolors
    }
    ]
    
    with open(namecmap+".json", "w") as outfile:
        json.dump(json_cmap, outfile)

json_rand_dictionary(Ncolors, 'RandomCmap', first_color_black=True)
