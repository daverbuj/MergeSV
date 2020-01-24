# MergeSV	MergeSV
```	
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

MergeSV is a command line tool used merge BED file tracks based on both reciprocal overlap and proximity to the ends of the tracks. 	

## Installation	

Use the package manager [pip](https://pip.pypa.io/en/stable/) to install mergesv.	

```bash	
pip install mergesv	
```	

## Usage	

```	
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
The merge is done by ensuring that the tracks have a reciprocal overlap (default 80%) as well as proximity to the starts and ends (default 1000 bp). The latter is particularly important in cases of large tracks where they share 80% overlap, but the start and end values differ largely. The merge is performed recursively, guaranteeing all tracks are adequately merged.	

The first merged track will be selected based on the median value of the group of tracks, but the start and end value of subsequent merges will be decided based on the previously-merged track with the largest number of merges.	

CIPOS and CIEND encompass the range of every track that becomes merged.	

## Contributing	
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.	

Please make sure to update tests as appropriate.	

## License	
[MIT](https://choosealicense.com/licenses/mit/)
