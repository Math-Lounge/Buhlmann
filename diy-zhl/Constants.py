
# ZH-L12 from "Decompression"
#
ZHL12N = {
    1 :  { "t" : 2.65,  "a" : 2.2,   "b" : 0.82 },
    2 :  { "t" : 7.94,  "a" : 1.5,   "b" : 0.82 },
    3 :  { "t" : 12.2,  "a" : 1.08,  "b" : 0.825 },
    4 :  { "t" : 18.5,  "a" : 0.9,   "b" : 0.835 },
    5 :  { "t" : 26.5,  "a" : 0.75,  "b" : 0.845 },
    6 :  { "t" : 37.0,  "a" : 0.58,  "b" : 0.86 },
    7 :  { "t" : 53.0,  "a" : 0.47,  "b" : 0.87 },
    8 :  { "t" : 79.0,  "a" : 0.45,  "b" : 0.89 },
    9 :  { "t" : 114.0, "a" : 0.45,  "b" : 0.89 },
    10 : { "t" : 146.0, "a" : 0.455, "b" : 0.934 },
    11 : { "t" : 185.0, "a" : 0.455, "b" : 0.934 },
    12 : { "t" : 238.0, "a" : 0.38,  "b" : 0.944 },
    13 : { "t" : 304.0, "a" : 0.255, "b" : 0.962 },
    14 : { "t" : 397.0, "a" : 0.255, "b" : 0.962 },
    15 : { "t" : 503.0, "a" : 0.255, "b" : 0.962 },
    16 : { "t" : 635.0, "a" : 0.255, "b" : 0.962 },
}

ZHL12He = {
    1 :  { "t" : 1.0,   "a" : 2.2,   "b" : 0.82 },
    2 :  { "t" : 3.0,   "a" : 1.5,   "b" : 0.82 },
    3 :  { "t" : 4.6,   "a" : 1.08,  "b" : 0.825 },
    4 :  { "t" : 7.0,   "a" : 0.9,   "b" : 0.835 },
    5 :  { "t" : 10.0,  "a" : 0.75,  "b" : 0.845 },
    6 :  { "t" : 14.0,  "a" : 0.58,  "b" : 0.86 },
    7 :  { "t" : 20.0,  "a" : 0.47,  "b" : 0.87 },
    8 :  { "t" : 30.0,  "a" : 0.45,  "b" : 0.89 },
    9 :  { "t" : 43.0,  "a" : 0.45,  "b" : 0.89 },
    10 : { "t" : 55.0,  "a" : 0.515, "b" : 0.926 },
    11 : { "t" : 70.0,  "a" : 0.515, "b" : 0.926 },
    12 : { "t" : 90.0,  "a" : 0.515, "b" : 0.926 },
    13 : { "t" : 115.0, "a" : 0.515, "b" : 0.926 },
    14 : { "t" : 150.0, "a" : 0.515, "b" : 0.926 },
    15 : { "t" : 190.0, "a" : 0.515, "b" : 0.926 },
    16 : { "t" : 240.0, "a" : 0.515, "b" : 0.926 },
}

# ZH_L16: several sources incl. a photo of a page from Tauchmedizin @
# http://www.nigelhewitt.co.uk/stuff/aab.jpg
# It appears nobody has Helium numbers for "-A" and "-C", nor for the 4-minute TC
# 5-minute TC is keyed as 1.1.
# Why does the 'a' column have multiple values? Big-X corresponds to ZHL-16_*X*
# We primarily care about the 'C' model. We can observe the differences in the image linked
ZHL16N_Base = {
    2   : { "t" : 8.0,   "b" : 0.6514, "a" : { "A" : 1.0,    "B" : 1.0,    "C" : 1.0 } },
    3   : { "t" : 12.5,  "b" : 0.7222, "a" : { "A" : 0.8618, "B" : 0.8618, "C" : 0.8618 } },
    4   : { "t" : 18.5,  "b" : 0.7825, "a" : { "A" : 0.7562, "B" : 0.7562, "C" : 0.7562 } },
    5   : { "t" : 27.0,  "b" : 0.8126, "a" : { "A" : 0.6667, "B" : 0.6667, "C" : 0.62 } },
    6   : { "t" : 38.3,  "b" : 0.8434, "a" : { "A" : 0.5933, "B" : 0.56,   "C" : 0.5043 } },
    7   : { "t" : 54.3,  "b" : 0.8693, "a" : { "A" : 0.5282, "B" : 0.4947, "C" : 0.441 } },
    8   : { "t" : 77.0,  "b" : 0.891,  "a" : { "A" : 0.4701, "B" : 0.45,   "C" : 0.4 } },
    9   : { "t" : 109.0, "b" : 0.9092, "a" : { "A" : 0.4187, "B" : 0.4187, "C" : 0.375 } },
    10  : { "t" : 146.0, "b" : 0.9222, "a" : { "A" : 0.3798, "B" : 0.3798, "C" : 0.35 } },
    11  : { "t" : 187.0, "b" : 0.9319, "a" : { "A" : 0.3497, "B" : 0.3497, "C" : 0.3295 } },
    12  : { "t" : 239.0, "b" : 0.9403, "a" : { "A" : 0.3223, "B" : 0.3223, "C" : 0.3065 } },
    13  : { "t" : 305.0, "b" : 0.9477, "a" : { "A" : 0.2971, "B" : 0.285,  "C" : 0.2835 } },
    14  : { "t" : 390.0, "b" : 0.9544, "a" : { "A" : 0.2737, "B" : 0.2737, "C" : 0.261 } },
    15  : { "t" : 498.0, "b" : 0.9602, "a" : { "A" : 0.2523, "B" : 0.2523, "C" : 0.248 } },
    16  : { "t" : 635.0, "b" : 0.9653, "a" : { "A" : 0.2327, "B" : 0.2327, "C" : 0.2327 } },
}

ZHL16N_4m = {
    1   : { "t" : 4.0,   "b" : 0.505,  "a" : { "A" : 1.2599, "B" : 1.2599, "C" : 1.2599 } },
    ** ZHL16N_Base,
}

ZHL16N_5m = {
    1.1 : { "t" : 5.0,   "b" : 0.5578, "a" : { "A" : 1.1696, "B" : 1.1696, "C" : 1.1696 } },
    ** ZHL16N_Base,
}

ZHL16N = {
    1   : { "t" : 4.0,   "b" : 0.505,  "a" : { "A" : 1.2599, "B" : 1.2599, "C" : 1.2599 } },
    1.1 : { "t" : 5.0,   "b" : 0.5578, "a" : { "A" : 1.1696, "B" : 1.1696, "C" : 1.1696 } },
    ** ZHL16N_Base,
}

# Note that "B" is a misnomer as some implementations call it "C",
# some call it "B", and nobody can read German and/or drop a hundred
# bucks on Tauchmedizin to see what Herr Buhlmann actually said.
#
ZHL16He = {
    1  :  { "t" : 1.51,   "a" : { "B" : 1.7424 }, "b" : 0.4245 },
    1.1 : { "t" : 1.88,   "a" : { "B" : 1.6189 }, "b" : 0.477 },
    2 :   { "t" : 3.02,   "a" : { "B" : 1.383 },  "b" : 0.5747 },
    3 :   { "t" : 4.72,   "a" : { "B" : 1.1919 }, "b" : 0.6527 },
    4 :   { "t" : 6.99,   "a" : { "B" : 1.0458 }, "b" : 0.7223 },
    5 :   { "t" : 10.21,  "a" : { "B" : 0.922 },  "b" : 0.7582 },
    6 :   { "t" : 14.48,  "a" : { "B" : 0.8205 }, "b" : 0.7957 },
    7 :   { "t" : 20.53,  "a" : { "B" : 0.7305 }, "b" : 0.8279 },
    8 :   { "t" : 29.11,  "a" : { "B" : 0.6502 }, "b" : 0.8553 },
    9 :   { "t" : 41.20,  "a" : { "B" : 0.595 },  "b" : 0.8757 },
    10 :  { "t" : 55.19,  "a" : { "B" : 0.5545 }, "b" : 0.8903 },
    11 :  { "t" : 70.69,  "a" : { "B" : 0.5333 }, "b" : 0.8997 },
    12 :  { "t" : 90.34,  "a" : { "B" : 0.5189 }, "b" : 0.9073 },
    13 :  { "t" : 115.29, "a" : { "B" : 0.5181 }, "b" : 0.9122 },
    14 :  { "t" : 147.42, "a" : { "B" : 0.5176 }, "b" : 0.9171 },
    15 :  { "t" : 188.24, "a" : { "B" : 0.5172 }, "b" : 0.9217 },
    16 :  { "t" : 240.03, "a" : { "B" : 0.5119 }, "b" : 0.9267 },
}

# From Deco for Divers
# Do not use: 
#  need to convert to bar and delta-M to use,
#  M0 is in msw and M = Delta M * Depth + M0,
#
WORKMAN = {
    1 : { "t" : 5.0, "M0" : 31.5, "M" : 1.8 },
    2 : { "t" : 10.0, "M0" : 26.8, "M" : 1.6 },
    3 : { "t" : 20.0, "M0" : 21.9, "M" : 1.5 },
    4 : { "t" : 40.0, "M0" : 17.0, "M" : 1.4 },
    5 : { "t" : 80.0, "M0" : 16.4, "M" : 1.3 },
    6 : { "t" : 120.0, "M0" : 15.8, "M" : 1.2 },
    7 : { "t" : 160.0, "M0" : 15.5, "M" : 1.15 },
    8 : { "t" : 200.0, "M0" : 15.5, "M" : 1.1 },
    9 : { "t" : 240.0, "M0" : 15.2, "M" : 1.1 },
}

# Also from Deco for Divers
#
# Since DSAT's primary concern is no-stop diving, it only uses M0
#  -- there is no Delta M i.e. Delta M = 1
# values are in msw
#
# Convert to Buhlmann with
#  m_w2b( M0 = NNN / 10, dM = 1, P = 1 )
# and run ZHL with DSAT compartments and M-values
#
DSAT = {
    1 : { "t" : 5.0, "M0" : 30.42 },
    2 : { "t" : 10.0, "M0" : 25.37 },
    3 : { "t" : 20.0, "M0" : 20.54 },
    4 : { "t" : 30.0, "M0" : 18.34 },
    5 : { "t" : 40.0, "M0" : 17.11 },
    6 : { "t" : 60.0, "M0" : 15.79 },
    7 : { "t" : 80.0, "M0" : 15.11 },
    8 : { "t" : 100.0, "M0" : 14.69 },
    9 : { "t" : 120.0, "M0" : 14.41 },
    10 : { "t" : 160.0, "M0" : 14.06 },
    11 : { "t" : 200.0, "M0" : 13.84 },
    12 : { "t" : 240.0, "M0" : 13.69 },
    13 : { "t" : 360.0, "M0" : 13.45 },
    14 : { "t" : 480.0, "M0" : 13.33 },
}