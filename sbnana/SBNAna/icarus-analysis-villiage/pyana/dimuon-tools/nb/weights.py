cv = [
  ("ppfx", "cv")
]

fluxsyst = [
  "beam_div",
  "beam_shift_x",
  "beam_spot",
  "horn1_x",
  "horn1_y",
  "horn_current_plus",
  "water_layer",
] + ["pca%i" % i for i in range(20)]

fluxsyst = [(f, s) for f in fluxsyst for s in ("ps", "ms")]
