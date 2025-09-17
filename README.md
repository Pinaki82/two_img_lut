# two_img_lut

Generate a 3D LUT (.cube) file by comparing two images: a source image and a colour-graded target image.

---

## Authorship & License

**License**: MIT-0 (Public Domain equivalent)  
**Written by**: ChatGPT on 17 September 2025  
**Copyright**: Pinaki Gupta

---

## Features

- Accepts two images: the original (source) and the colour-graded version (target).  
- Produces a `.cube` 3D LUT mapping source → target.  
- Supports a default LUT size if not specified.  
- Prints elapsed wall-clock time after processing.  
- Linux only (uses `clock_gettime`, etc.).

---

## Requirements

- Linux operating system.  
- GCC (or a C compiler) supporting C11.  
- Math library (`libm`) available.  
- `stb_image.h` and `stb_image_write.h` [from [GitHub - nothings/stb: stb single-file public domain libraries for C/C++](https://github.com/nothings/stb.git)] placed under `external/stb/` (or wherever included via the Makefile).

---

## Installation

Assuming your directory structure is:

```css
.
your-project/
├── external/
│ └── stb/
│ ├── stb_image.h
│ └── stb_image_write.h
├── src/
│ └── two_img_lut.c
│ └── lutgen.h
├── Makefile
└── README.md
```

To build and install:

Clone the repository (if you haven’t already). You can use the --recursive flag so submodules are initialised at clone time:

```bash
git clone --recursive https://github.com/Pinaki82/two_img_lut.git
cd two_img_lut/
```

If you already have the repo (without submodules), then after pulling, run:

```bash
git submodule update --init --recursive
```

After that, to update the main repo and submodules later, you might want to do:

```bash
git pull --recurse-submodules
git submodule update --init --recursive
```

Or, configure git so that git pull always recurses into submodules.

Then,

```bash
make
make install
```

This compiles the `two_img_lut binary` and copies it to `~/.local/bin`.
Make sure `~/.local/bin` is in your `PATH`:

```bash
export PATH="$HOME/.local/bin:$PATH"
```

Put that into your shell startup (e.g. `~/.bashrc` or `~/.zshrc`) if not already.

```bash
two_img_lut <source.png> <target.png> [lut_size[optional]] <output.cube>
```

- If lut_size is omitted, a default size (e.g. 33) is used.

- Examples:

```bash
$ ls
source.png
source_corrected.png

# Generate LUT with default size
two_img_lut source.png source_corrected.png corrected_lut.cube

# Or specify LUT size
two_img_lut source.png source_corrected.png 17 small_lut.cube
```

Input images must have the same dimensions.

**Notes & Tips**

- The generated .cube file is ASCII, assumes RGB values in the [0,1] domain for both input & output.

- Cells in the LUT with no sample data are filled via nearest-neighbour fallback.

- A higher lut_size gives finer colour mapping but slower processing and a larger LUT file.

- Ensure both source and target images are in the same colour space/gamma curve before generating the LUT; mismatched gamma can lead to incorrect results.

```css
.
├── external/
│   └── stb/
│       ├── stb_image.h
│       └── stb_image_write.h
├── src/
│   └── two_img_lut.c
│   └── lutgen.h
├── Makefile
└── README.md
```

**License**

MIT-0 — Public Domain equivalent. See source file for full text.

```css
MIT No Attribution

Copyright <year> <copyright holders>

Permission is hereby granted, free of charge, to any person obtaining a copy of this
software and associated documentation files (the "Software"), to deal in the Software
without restriction, including without limitation the rights to use, copy, modify,
merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
```
