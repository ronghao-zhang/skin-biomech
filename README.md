# skin-biomech-scripts
This repo contains MatLab scripts that aims to process skin biomechanics signals &amp; to reveal the neural basis for somatosensation 

## Code & Scripts
### The directory `fabricScan`
This directory uses height maps of fabric surfaces (Scanned in 2021) as inputs. The `fs1_initialziation` script defines parameters that will be used later. 
The `fs2_preprocessing` script perfrom preprocessing on the height maps. The preprocess include downsampling, detrending, etc. 
The `fs3_preindentation+afferent_response` script simulates the response of afferents and calculate a customized pre-indentation value for each texture since we want to ensure that the force applied on each should be the same.   
Once the customized pre-indentation is obtained, change the corresponding part in the `fs2_preprocessing` part and re-run `fs2_*` and `fs3_*`.

## The directory `freqAnalysis`
This directory contains scripts that can calibrate the height of, preprocess, and plot height maps. The discription of each script can be found in the first few lines of each scripts. 

## The directory `others`
This directory contains scripts than can automatically create height map of square gratings. 

## Results & Figures

1. **Figure_01**: Downsampled Height Maps of Fabric Surface Patterns using a Low Resolution of 0.1 mm/pixel.
2. **Figure_02**: Detrended & Downsampled Height Maps of Fabric Surface Patterns using a Low Resolution of 0.1 mm/pixel.
3. **Figure_03**: The Force Applied Over Time using an universal pre-indentation for textures. Blue line is the force at each time point, Red Line is the average force over the period.
4. **Figure_04**: Detrended & Downsampled Height Maps of Fabric Surface Patterns after applied a customized pre-indentation value to make sure that the pressure is 0.05 Newtons per mm of indentation with a 2mm diameter pin. 
5. **Figure_05**: The spike train plot that shows the firing over time of different afferent neurons when sliding fingertip on the texture using customized pre-indentation.
6. **Figure_06**: The Force Applied Over Time after appying a customized pre-indentation (just a check). 
7. **Figure_07**: The Square Gratings created for testing the skin indentation. 
8. **Figure_08**: The linear relationship between the height of the gel and the height of the original pattern. 
9. **Results_of_Cross_Correlation**: All matched pairs of gelsight and ground truth. 
