files_my_gen = [string for string in open('noProbeFilterDecayFilter_MiniAOD.txt').readlines() if string != ""]

for f in files_my_gen:
    print f