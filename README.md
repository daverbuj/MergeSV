# MergeSV
```bash
o     o                            .oPYo. o     o 
8b   d8                            8      8     8 
8`b d'8 .oPYo. oPYo. .oPYo. .oPYo. `Yooo. 8     8 
8 `o' 8 8oooo8 8  `' 8    8 8oooo8     `8 `b   d' 
8     8 8.     8     8    8 8.          8  `b d'  
8     8 `Yooo' 8     `YooP8 `Yooo' `YooP'   `8'   
<------------>            8     <----------->     
  <------------>       ooP'  <----------------->  
      <--------->      <-------->        <------->
```

MergeSV is a command line tool 

## Installation

Use the package manager [pip](https://pip.pypa.io/en/stable/) to install mergesv.

```bash
pip install mergesv
```

## Usage

```bash
mergesv <-i BED> [options]

Required Arguments:
	-i    FILE    BED file

Options:
	-o    FILE    output file [bed-file.mergesv.bed]
	-r    FLOAT   reciprocal overlap [0.8]
	-s    INT     distance (bp) from start and end values to be considered a merge [1000]
	-c            Include number of merged SVs in output
```

## How the merge works
This is how the merge works

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License
[MIT](https://choosealicense.com/licenses/mit/)
