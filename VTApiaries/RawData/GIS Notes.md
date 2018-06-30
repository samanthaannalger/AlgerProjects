GIS Notes:



To display X Y coordinates:



Data in Excel doc

File -> Add Data -> Add XY Data 

​	Select the excel Doc, x field is the long and y field is the lat, Change projected coordinates to WGS 1984



**To create a thematic map (counties color coded by hive losses, for example)**

***There are problems if the excel doc has 'NA's in the attribute column you are using for the thematic map: i.e.: colony loss, Solution: strip NA's from excel doc first before adding XY data to ArcGIS.

1. Make the xy layer into a shape file, data->export-> rename as shape file

2. Need: county Boundary Polygon layer

3. Use the shape file of x y data to create a raster data layer: Feature to Raster tool: 
   Input features: the shape file
   Field: colony loss column from the xy data set
   Output: name the layer

   ​		Results: is a raster layer that can be used in Zonal Statistics...


4. Search for the tool: Zonal Statistics as a table:
   Input raster or feature zone: County polygons (zone field will be object ID)Input raster: Raster layer created above
   name output table
   statistics type: ALL
5. Join the polygon layer to the Zonal Statistics table 
   Right Click county boundary polygon label -> Joins and relates->Join-> ZonalStatsTable


6. After joined, change symbology so that each county is colored accordingly….
   Go to properties of county boundary polygon label-> Symbology-> Categories 
   change 'value field' to Colony Loss, unselect "all other values", click, "add all values", Change color ramp accordingly