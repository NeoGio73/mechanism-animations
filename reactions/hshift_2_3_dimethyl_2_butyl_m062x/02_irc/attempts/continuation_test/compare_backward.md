# Continuation test, hshift_2_3_dimethyl_2_butyl_m062x, backward side

Continuation from the last frame of `attempts/open_side_stopped_short_1218/backward` vs rerun from the TS (`true_irc/backward`), compared at equal arc length s (rerun interpolated). Criteria: each new frame heavy-atom RMSD ≤ 0.05 Å and |ΔE| ≤ 0.1 kJ/mol; both ends within 0.5 kJ/mol of hshift_tertiary/open_cation, end to end |ΔE| ≤ 0.1 kJ/mol and heavy-atom RMSD ≤ 0.2 Å. Early comparison (J. Picione, 2026-11-16): both runs stopped once the rerun had 50 frames past the old side's last frame; verdict = the per-frame criterion over at least 50 compared new frames; the end criterion is not applied. **Agree: yes**

- Earlier frames (112, same settings from the TS): max heavy-atom RMSD 0.0000 Å, max |ΔE| 0.0000 kJ/mol.
- Start frame recomputed: +0.0000 kJ/mol.
- New frames compared: 54 (beyond the rerun's end: 98); heavy-atom RMSD max 0.0000 Å, median 0.0000 Å; max |ΔE| 0.0001 kJ/mol.
- Ends: rerun 167 frames, s = 8.34, 3.136 kJ/mol above the minimum; continuation 264 frames (152 new), s = 13.29, 1.468 kJ/mol above; end to end heavy-atom RMSD 0.236 Å, ΔE -1.668 kJ/mol.

| frame | s (bohr·amu½) | heavy-atom RMSD (Å) | ΔE continuation − rerun (kJ/mol) | E − minimum (kJ/mol) |
|---|---|---|---|---|
| 113 | 5.587 | 0.0000 | +0.0000 | 4.106 |
| 114 | 5.637 | 0.0000 | -0.0000 | 4.084 |
| 115 | 5.688 | 0.0000 | -0.0000 | 4.064 |
| 116 | 5.738 | 0.0000 | -0.0000 | 4.043 |
| 117 | 5.788 | 0.0000 | -0.0000 | 4.023 |
| 118 | 5.839 | 0.0000 | -0.0000 | 4.003 |
| 119 | 5.889 | 0.0000 | -0.0000 | 3.983 |
| 120 | 5.940 | 0.0000 | -0.0000 | 3.963 |
| 121 | 5.991 | 0.0000 | -0.0000 | 3.944 |
| 122 | 6.041 | 0.0000 | -0.0000 | 3.924 |
| 123 | 6.092 | 0.0000 | -0.0000 | 3.904 |
| 124 | 6.142 | 0.0000 | -0.0000 | 3.884 |
| 125 | 6.193 | 0.0000 | -0.0000 | 3.863 |
| 126 | 6.243 | 0.0000 | -0.0000 | 3.843 |
| 127 | 6.294 | 0.0000 | +0.0000 | 3.822 |
| 128 | 6.345 | 0.0000 | +0.0000 | 3.801 |
| 129 | 6.395 | 0.0000 | +0.0000 | 3.780 |
| 130 | 6.446 | 0.0000 | -0.0000 | 3.758 |
| 131 | 6.496 | 0.0000 | +0.0000 | 3.736 |
| 132 | 6.547 | 0.0000 | -0.0000 | 3.714 |
| 133 | 6.598 | 0.0000 | -0.0000 | 3.692 |
| 134 | 6.648 | 0.0000 | -0.0000 | 3.670 |
| 135 | 6.699 | 0.0000 | -0.0000 | 3.648 |
| 136 | 6.750 | 0.0000 | -0.0000 | 3.625 |
| 137 | 6.801 | 0.0000 | -0.0000 | 3.603 |
| 138 | 6.851 | 0.0000 | -0.0000 | 3.581 |
| 139 | 6.902 | 0.0000 | -0.0000 | 3.559 |
| 140 | 6.953 | 0.0000 | -0.0000 | 3.537 |
| 141 | 7.004 | 0.0000 | -0.0000 | 3.515 |
| 142 | 7.055 | 0.0000 | -0.0000 | 3.494 |
| 143 | 7.106 | 0.0000 | -0.0000 | 3.473 |
| 144 | 7.157 | 0.0000 | -0.0000 | 3.452 |
| 145 | 7.209 | 0.0000 | -0.0000 | 3.432 |
| 146 | 7.260 | 0.0000 | -0.0000 | 3.412 |
| 147 | 7.311 | 0.0000 | -0.0000 | 3.392 |
| 148 | 7.363 | 0.0000 | -0.0000 | 3.373 |
| 149 | 7.414 | 0.0000 | -0.0000 | 3.354 |
| 150 | 7.465 | 0.0000 | -0.0000 | 3.336 |
| 151 | 7.517 | 0.0000 | -0.0000 | 3.319 |
| 152 | 7.569 | 0.0000 | -0.0001 | 3.302 |
| 153 | 7.620 | 0.0000 | -0.0000 | 3.286 |
| 154 | 7.672 | 0.0000 | -0.0001 | 3.271 |
| 155 | 7.724 | 0.0000 | -0.0001 | 3.256 |
| 156 | 7.776 | 0.0000 | -0.0001 | 3.242 |
| 157 | 7.828 | 0.0000 | -0.0000 | 3.229 |
| 158 | 7.880 | 0.0000 | -0.0001 | 3.217 |
| 159 | 7.933 | 0.0000 | -0.0001 | 3.205 |
| 160 | 7.985 | 0.0000 | -0.0000 | 3.194 |
| 161 | 8.038 | 0.0000 | -0.0001 | 3.183 |
| 162 | 8.090 | 0.0000 | -0.0001 | 3.174 |
| 163 | 8.143 | 0.0000 | -0.0001 | 3.165 |
| 164 | 8.195 | 0.0000 | -0.0000 | 3.156 |
| 165 | 8.248 | 0.0000 | -0.0000 | 3.148 |
| 166 | 8.301 | 0.0000 | -0.0000 | 3.141 |
| 167 | 8.354 | beyond the rerun's end | — | — |
| 168 | 8.408 | beyond the rerun's end | — | — |
| 169 | 8.462 | beyond the rerun's end | — | — |
| 170 | 8.517 | beyond the rerun's end | — | — |
| 171 | 8.571 | beyond the rerun's end | — | — |
| 172 | 8.626 | beyond the rerun's end | — | — |
| 173 | 8.680 | beyond the rerun's end | — | — |
| 174 | 8.734 | beyond the rerun's end | — | — |
| 175 | 8.789 | beyond the rerun's end | — | — |
| 176 | 8.843 | beyond the rerun's end | — | — |
| 177 | 8.898 | beyond the rerun's end | — | — |
| 178 | 8.953 | beyond the rerun's end | — | — |
| 179 | 9.008 | beyond the rerun's end | — | — |
| 180 | 9.062 | beyond the rerun's end | — | — |
| 181 | 9.115 | beyond the rerun's end | — | — |
| 182 | 9.169 | beyond the rerun's end | — | — |
| 183 | 9.222 | beyond the rerun's end | — | — |
| 184 | 9.275 | beyond the rerun's end | — | — |
| 185 | 9.327 | beyond the rerun's end | — | — |
| 186 | 9.380 | beyond the rerun's end | — | — |
| 187 | 9.432 | beyond the rerun's end | — | — |
| 188 | 9.484 | beyond the rerun's end | — | — |
| 189 | 9.536 | beyond the rerun's end | — | — |
| 190 | 9.587 | beyond the rerun's end | — | — |
| 191 | 9.639 | beyond the rerun's end | — | — |
| 192 | 9.688 | beyond the rerun's end | — | — |
| 193 | 9.738 | beyond the rerun's end | — | — |
| 194 | 9.788 | beyond the rerun's end | — | — |
| 195 | 9.838 | beyond the rerun's end | — | — |
| 196 | 9.888 | beyond the rerun's end | — | — |
| 197 | 9.939 | beyond the rerun's end | — | — |
| 198 | 9.989 | beyond the rerun's end | — | — |
| 199 | 10.040 | beyond the rerun's end | — | — |
| 200 | 10.091 | beyond the rerun's end | — | — |
| 201 | 10.142 | beyond the rerun's end | — | — |
| 202 | 10.193 | beyond the rerun's end | — | — |
| 203 | 10.243 | beyond the rerun's end | — | — |
| 204 | 10.294 | beyond the rerun's end | — | — |
| 205 | 10.344 | beyond the rerun's end | — | — |
| 206 | 10.395 | beyond the rerun's end | — | — |
| 207 | 10.445 | beyond the rerun's end | — | — |
| 208 | 10.496 | beyond the rerun's end | — | — |
| 209 | 10.546 | beyond the rerun's end | — | — |
| 210 | 10.596 | beyond the rerun's end | — | — |
| 211 | 10.646 | beyond the rerun's end | — | — |
| 212 | 10.696 | beyond the rerun's end | — | — |
| 213 | 10.746 | beyond the rerun's end | — | — |
| 214 | 10.796 | beyond the rerun's end | — | — |
| 215 | 10.845 | beyond the rerun's end | — | — |
| 216 | 10.895 | beyond the rerun's end | — | — |
| 217 | 10.945 | beyond the rerun's end | — | — |
| 218 | 10.995 | beyond the rerun's end | — | — |
| 219 | 11.044 | beyond the rerun's end | — | — |
| 220 | 11.094 | beyond the rerun's end | — | — |
| 221 | 11.143 | beyond the rerun's end | — | — |
| 222 | 11.192 | beyond the rerun's end | — | — |
| 223 | 11.242 | beyond the rerun's end | — | — |
| 224 | 11.291 | beyond the rerun's end | — | — |
| 225 | 11.340 | beyond the rerun's end | — | — |
| 226 | 11.390 | beyond the rerun's end | — | — |
| 227 | 11.439 | beyond the rerun's end | — | — |
| 228 | 11.489 | beyond the rerun's end | — | — |
| 229 | 11.538 | beyond the rerun's end | — | — |
| 230 | 11.588 | beyond the rerun's end | — | — |
| 231 | 11.637 | beyond the rerun's end | — | — |
| 232 | 11.687 | beyond the rerun's end | — | — |
| 233 | 11.737 | beyond the rerun's end | — | — |
| 234 | 11.786 | beyond the rerun's end | — | — |
| 235 | 11.836 | beyond the rerun's end | — | — |
| 236 | 11.886 | beyond the rerun's end | — | — |
| 237 | 11.936 | beyond the rerun's end | — | — |
| 238 | 11.986 | beyond the rerun's end | — | — |
| 239 | 12.036 | beyond the rerun's end | — | — |
| 240 | 12.086 | beyond the rerun's end | — | — |
| 241 | 12.136 | beyond the rerun's end | — | — |
| 242 | 12.186 | beyond the rerun's end | — | — |
| 243 | 12.236 | beyond the rerun's end | — | — |
| 244 | 12.286 | beyond the rerun's end | — | — |
| 245 | 12.336 | beyond the rerun's end | — | — |
| 246 | 12.386 | beyond the rerun's end | — | — |
| 247 | 12.436 | beyond the rerun's end | — | — |
| 248 | 12.486 | beyond the rerun's end | — | — |
| 249 | 12.536 | beyond the rerun's end | — | — |
| 250 | 12.586 | beyond the rerun's end | — | — |
| 251 | 12.636 | beyond the rerun's end | — | — |
| 252 | 12.686 | beyond the rerun's end | — | — |
| 253 | 12.736 | beyond the rerun's end | — | — |
| 254 | 12.786 | beyond the rerun's end | — | — |
| 255 | 12.837 | beyond the rerun's end | — | — |
| 256 | 12.887 | beyond the rerun's end | — | — |
| 257 | 12.937 | beyond the rerun's end | — | — |
| 258 | 12.987 | beyond the rerun's end | — | — |
| 259 | 13.037 | beyond the rerun's end | — | — |
| 260 | 13.087 | beyond the rerun's end | — | — |
| 261 | 13.137 | beyond the rerun's end | — | — |
| 262 | 13.187 | beyond the rerun's end | — | — |
| 263 | 13.237 | beyond the rerun's end | — | — |
| 264 | 13.287 | beyond the rerun's end | — | — |
