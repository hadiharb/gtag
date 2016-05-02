# gtag
gtag is an audio identification software. It lets you tag audio files and then identify unknown audio signals by comparing them to a known dataset. The technique is robust to noise such as loudspeaker - microphone audio capture.

Usage: gtag extract filename_audio filename_fingerprint
To extract the fingerprint of an audio file. It is very important to extract fingerprints from WAV, PCM, 8KHz, 16bits, mono audio files. The fingerprint filename must have the .fin extension.

Usage: gtag compare filename_audio directory_fingerprints
To compare an audio file (WAV, PCM, 8KHz, 16bits, mono) with the fingerprints (.fin files) in a directory.

Main contributor: Hadi Harb (hadi[DOT]harb[AT]ghanni[DOT]com) with the contribution of: Aliaksandr Paradzinets (aliaksandr[DOT]paradzinets[AT]ghanni[DOT]com)
