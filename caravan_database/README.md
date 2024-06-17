## Test of P model on caravn database
This read me is temporaty, I summed what I done so far
The caravan database is described in [Caravan paper](https://www.nature.com/articles/s41597-023-01975-w)
The data are avaible on [Zenodo](https://doi.org/10.5281/zenodo.7540792)

The script function_check.R search for the closest site of caravan database in respect to the tower site located on US.
I filter out the sites that are further of 0.05 degree and check if the total precipitation is the same (it should be since are closer)

In caravan database the streamflow is observed so I calculate the difference between precipitation and streamflow to get the AET.
I compare this AET with the AET obtained from P model
