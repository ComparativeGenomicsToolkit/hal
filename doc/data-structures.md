# HAL data structures

Note: this is currently a work in progress.  Right now, it is used to collect various
pieces of information that need to be explained in the full documentation.


## TopSegments

Each genome has a segment array that maps contiguous ranges of the genome to
ranges in the parent. The array is sorted by

### TopSegment entry
- startPosition
- parentIndex -
- reversed -
- bottomParseIndex - the index of the bottom segment array in the genome that contains the start coordinate of this top segment.
- bottomParseOffset
- paralogyIndex

### BottomSegments

### BottomSegment entry
- startPosition - 
- topParseIndex - 
