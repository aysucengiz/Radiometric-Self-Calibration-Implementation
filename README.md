# Radiometric-Self-Calibration-Implementation (Mitsunaga and Nayar 1999)

> Blog under the `blog` folder.
> 
### Dependencies
```pip install numpy matplotlib pillow```

### Running
```python main {FOLDER_PATH}```

Options: 
* `--custom`: Is the image set from the given examples or not. If not, a csv with the F values is expected with the same order as the images according to their names' alphabetical order.
* `--csv {EXPOSURE_CSV}`: Gives the path to said csv.

Example1:
```python main "./images/exposure sequences/exposure sequences/Nikon_Coolpix_E5400/" ```

Example2:
```python main "./images/mine/" --custom --csv="./images/mine/exposures.csv" ```

