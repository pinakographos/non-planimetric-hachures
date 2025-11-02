# Non-Planimetric Hachures
A QGIS & Blender-based method to generate automated hachure drawings of terrain. Like this:

<img width="3398" height="1230" alt="image" src="https://github.com/user-attachments/assets/d9401798-170a-4884-9923-e896d6eb0556" />


# Preamble
This method is a compantion to, and based off of, my planimetric [hachure method](https://github.com/pinakographos/Hachures).
This is version 1.0, which means that there are likely to be bugs; be warned. The two raster images in this repo (Uphill.png and Shading.png) are test images that have previously been known to work properly with the script (in QGIS 3.40.5-Bratislava), and should serve as a confirmation that you are using it correctly.

Before moving on, I want to offer my thanks to [Andrew Tyrrell](https://southarrowmaps.co.nz/), who made significant intellectual contributions to the method that is presented here, and helped with testing.

In order to use this script, you will need to first generate some images in Blender. And you will need a general understanding of what is going on. Much as it pains me to say, you have some homework to do. You will need to read this [draft article](https://docs.google.com/document/d/1GFSgxLI5F2UOevCXBG22GBKVbpuKddTl7MoBTrPHB3s/edit?usp=sharing) in order to properly implement the method.

# How to Run
This is a QGIS script, not a plugin, so there's no fancy user interface (I tried, but it turns out that making a plugin was going to be more than I could handle). But it's not too hard to get it running.
+ On this page, go up to the green "Code" button, click it, and select "Download ZIP".
+ Unpack the .ZIP file somewhere on your computer. It contains the script, as well as test images.
+ In QGIS, load in the images you want to work with (such as the test images).
+ Then, go to the "Plugins" menu and choose "Python Console." A window will pop up somewhere in your interface.
+ In the Python Console, there's a button that looks like a pencil and paper with the tooltip "Show Editor." Click that. <img width="63" alt="image" src="https://github.com/user-attachments/assets/17a99a69-331d-4d48-b441-f9371855f987" />
+ A new widow will open within the Python Console. In this window is a button that looks like a yellow folder with the tooltip "Open Script…." Click that. <img width="58" alt="image" src="https://github.com/user-attachments/assets/36b66c52-6c9f-49f4-8ca5-d9d7b5c8d24f" />
+ Browse to where you stored the script, and open it.
+ Now you can see the script contents. You can type in values for the user parameters (see below for a discussion).
+ When ready, press the "Run Script" button, which looks like a green arrow <img width="69" alt="image" src="https://github.com/user-attachments/assets/9327f315-7e40-47b6-ba6e-51ab7d654e6a" />. Note that there are **two** green arrows. You want the one that's in the window with the script.
+ Wait patiently for hachures to generate.

# Initial Parameters
Once you've prepared your Blender layers, you can load them into QGIS and then fill them in as the `shading_layer_name` and `upslope_layer_name` parameters.

The script comes with some default parameters, but you may choose to adjust them:
+ `min_hachure_spacing` and `max_hachure_spacing`: These specify how close or how far apart we'd like our hachures to be. The units are the pixel size of the DEM.
+ `stubFilter` cleans out tiny lines, and is set to `max_hachure_spacing`, but you can adjust the filter to be larger or smaller as you like.
+ `breaks_layer_name` lets you specify a vector layer of lines to avoid drawing hachures near.

# Advice
First off, **be patient**. This script can take a while to run, depending on the settings. While a 1000 × 1000px raster with a handful of hachures may process in seconds, if you want a detailed set of lines on a large terrain image, it could potentially run for a long time. Start small, and then work your way up to more detail and larger terrains once you get a sense of how long it will take.

Getting good results takes iteration. Try different settings, and try smoothing your DEM in different ways, and seeing how every adjustment affects the results. A manual hachure artist has a lot of decisions to make about where and how to draw lines. This script will do the drawing for you, but you still need to make the decisions about what parameters look best.

While I am happy to try and help troubleshoot, I do ask that you invest the time to first read through the [draft article](https://docs.google.com/document/d/1GFSgxLI5F2UOevCXBG22GBKVbpuKddTl7MoBTrPHB3s/edit?usp=sharing) to ensure that you understand the methdology, before filing a bug report or asking questions. In truth, I am rather a novice programmer, so my ability to aid you is limited.

# Still Stuck?
If you're interested in achieving this hachure style, but don't want to wade through the documentation, code, etc: you can also [hire me](mailto:somethingaboutmaps@gmail.com) to make a terrain image for you to use in your projects!
