# Smooth scalars for VTK mesh 

This program smooths the point data scalars based on the maximum distance between two extreme points on the mesh 

Note that this maximum distance ``hard-coded`` is used to calculate the weight or influence of each point. The weight is based on the proportion of the distance between points.  

It is strongly recommended that this distance value is adapted for each mesh. 

## Usage 
The usage is through command line: 
```
SmoothScalars <input_VTK> <output_VTK> <iterations>
```

## Parameters

The first two parameters ```input_VTK``` and ```output_VTK``` are the input and output VTK files. The number of smoothing ```iterations``` is dependent on the application 


## Author 
Note that the code was adapted from [this VTK users forum post](http://vtk.1045678.n5.nabble.com/How-to-smooth-or-filter-polydata-scalars-td1243842.html) and the credits go to Fabio Meneghini. 

```
Dr. Rashed Karim 
Department of Biomedical Engineering 
King's College London 
```
