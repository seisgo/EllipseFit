## Workflow of ellipse fitting through matlab files
### 1. Prepare examples of scatter data points to perform ellipse fitting
>- Matlab file: prepareEllipseData.m  
>- Description: 4 types of test data:
 (1) hyperbola shaped scatter data; 
 (2) hyperbola or ellipse shaped scatter data;
 (3) ideal ellipse shaped scatter data;
 (4) noisy ellipse shaped scatter data;
>- Output: (1) hyperbolaData.mat; (2) hyperEllipData.mat;
           (3) ellipseData.mat;   (4) noisyEllipData.mat;

### 2. Do RANSAC filtering before ellipse fitting for noisy data (_Optional_)
>- Matlab file: ellipseDataFilter_RANSAC.m
>- Type: function
>- Usage: called by main program

### 3. Do ellipse fitting 
>- Description: 5 types of functions provided
>- Note: some following methods may fail in certain circumstances
>- Type: functions
>- **funcEllipseFit_nlinfit.m**: _nlinfit_ for conic section
>- **funcEllipseFit_OGal.m**: function provided by _Ohad Gal_
>- **funcEllipseFit_BFisher.m**: direct fit function provided by _Bob Fisher_
>- **funcEllipseFit_direct.m**: direct fit function (same as Bob Fisher's) provided by _Nikolai Chernov_
>- **funcEllipseFit_RBrown.m**: function provided by _Richard Brown_

### 4. Plot fitted ellipse (_For certain case_)
>- Matlab file: plotellipse.m
>- Type: function
>- Usage: called by main program
>- Note: particularly for fitted ellipse by _funcEllipseFit_RBrown_

### 5. Main program calling above functions and utilizing above data
>- Matlab file: doEllipseFit.m
>- Description: (1) load saved data; (2) choose whether do RANSAC filtering;
                (3) do ellipse fitting; (4) plot fitted ellipse and make some comparison

### Appendix
>- Matlab file: compareEllipsePlotForm.m
>- Description: compare ellipse plot form as plot and plotm, 
                where 'plot' plots ellipse in Cartesian coordinate;
                while 'plotm' plots ellipse in map coordinate.
>- Html file: compareEllipsePlotForm.html, compareEllipsePlotForm.png 
              and compareEllipsePlotForm_01.png in html directory.