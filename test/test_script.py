files_my_gen = [string for string in open('noProbeFilterDecayFilter_MiniAOD.txt').readlines() if len(string) > 10]

for f in files_my_gen:
    if len(f) < 10:
        print f

