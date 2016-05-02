all: gtag_open
gtag_open: gtag_open.cpp GFingerPrint.h GFingerPrint.cpp GFPExtractor.h GFPExtractor.cpp GFFT.cpp GFFT.h
	g++ gtag_open.cpp GFingerPrint.cpp GFPExtractor.cpp GFFT.cpp -o gtag
