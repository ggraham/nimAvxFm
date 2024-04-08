import cligen, times

{.push header: "<AwFmIndex.h>".}
type 
  AwFmAlphabetType = enum
    AwFmAlphabetAmino = 1,
    AwFmAlphabetDna = 2,
    AwFmAlphabetRna = 3
  AwFmIndexConfiguration {.importc: "struct AwFmIndexConfiguration".} = object
    suffixArrayCompressionRatio: uint8
    kmerLengthInSeedTable: uint8
    alphabetType: AwFmAlphabetType
    keepSuffixArrayInMemory: bool
    storeOriginalSequence: bool
  AwFmIndex {.importc: "struct AwFmIndex".} = object
  AwFmSearchRange {.importc: "struct AwFmSearchRange".} = object
  AwFmReturnCode = enum
    AwFmFileAlreadyExists = -15, AwFmErrorSuffixArrayNull = -14,
    AwFmErrorDbSequenceNull  = -13, AwFmFileWriteFail = -12,
    AwFmFileReadFail = -11, AwFmFileOpenFail = -10,
    AwFmFileFormatError = -9, AwFmNoDatabaseSequenceGiven = -8,
    AwFmNoFileSrcGiven = -7, AwFmIllegalPositionError = -6,
    AwFmSuffixArrayCreationFailure = -5, AwFmNullPtrError = -4,
    AwFmAllocationFailure = -3, AwFmUnsupportedVersionError = -2,
    AwFmGeneralFailure = -1, AwFmSuccess = 1, AwFmFileReadOkay = 2,
    AwFmFileWriteOkay = 3
  AwFmKmerSearchData {.importc: "struct AwFmKmerSearchData".} = object
    kmerString: cstring
    kmerLength: uint64
    positionList: uint64
    count: uint32
    capacity: uint32
  AwFmKmerSearchList {.importc: "struct AwFmKmerSearchList".} = object
    capacity: csize_t
    count: csize_t
    kmerSearchData: array[100000001, AwFmKmerSearchData]

proc awFmCreateIndexFromFasta(index: ptr ref AwFmIndex,
                              config: ptr AwFmIndexConfiguration,
                              fastaSrc: cstring,
                              indexFileSrc: cstring): AwFmReturnCode {.importc.}

proc awFmReadIndexFromFile(index: ptr ref AwFmIndex,
                           file: cstring,
                           keepSuffixArrayInMemory: bool = false): AwFmReturnCode {.importc.}

proc awFmCreateKmerSearchList(capacity: csize_t): ref AwFmKmerSearchList {.importc.}

proc awFmDeallocIndex(index: ref AwFmIndex) {.importc.}
proc awFmDeallocKmerSearchList(searchList: ref AwFmKmerSearchList) {.importc.}

proc awFmFindSearchRangeForString(index: ptr AwFmIndex,
                                  kmer: cstring,
                                  kmerLength: uint16): AwFmSearchRange {.importc.}
proc awFmGetNumSequences(index: ref AwFmIndex): cuint {.importc.}
proc awFmSingleKmerExists(index: ref AwFmIndex, km: cstring, kmerLength: uint16): bool {.importc.}
proc awFmParallelSearchCount(index: ref AwFmIndex, searchList: ref AwFmKmerSearchList, numThreads: uint8) {.importc.}
{.pop.}

proc map(index: string, km: string = "AGCCTGA", reps: int = 1, threads: int = 4) =
  let
    indexP = new(AwFmIndex)
  echo awFmReadIndexFromFile(indexP.unsafeAddr, index, true)

  var
    searchList = awFmCreateKmerSearchList(reps.csize_t)
  for i in 0..reps:
    searchList.kmerSearchData[i].kmerString = km
    searchList.kmerSearchData[i].kmerLength = km.len.uint64
  searchList.count = reps.csize_t
  let time = cpuTime()
  awFmParallelSearchCount(indexP,
                          searchList,
                          4)
  echo "Time parallelSearchCount: ", cpuTime() - time
  echo searchList.kmerSearchData[0].count
  awFmDeallocIndex(indexP)

proc index(fasta: string,
           index: string,
           saCompressionRatio: uint8 = 2,
           kmerLength: uint8 = 12) =
  let 
    indexP = new(AwFmIndex)
    conf = AwFmIndexConfiguration(
      suffixArrayCompressionRatio: saCompressionRatio,
      kmerLengthInSeedTable: kmerLength,
      alphabetType: AwFmAlphabetDna,
      keepSuffixArrayInMemory: false,
      storeOriginalSequence: false)
  echo awFmCreateIndexFromFasta(indexP.unsafeAddr, conf.unsafeAddr, fasta, index)
  awFmDeallocIndex(indexP)

dispatchMulti([map], [index])

