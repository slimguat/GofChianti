;===============================================================
; File: chi_find_transition.pro
; Find the closest line (by wavelength) in a CHIANTI .wgfa file
; Usage:
;   r = chi_find_transition('o_5', 764.0, /VERBOSE)
;   r = chi_find_transition('ne__', 770.0, CHIANTI_ROOT='/path/to/xuvtop', /VERBOSE)
;===============================================================

; Make "./tmp/ioneq_name_YYYYMMDDHHMMSS/" in pure IDL
FUNCTION datetime_stamp
  COMPILE_OPT idl2
  jd = SYSTIME(/JULIAN)
  CALDAT, jd, month, day, year, hour, minute, second

  ; Format: YYYYMMDDHHMMSS
  stamp = STRING(year, FORMAT='(I4.4)') + $
          STRING(month, FORMAT='(I2.2)') + $
          STRING(day, FORMAT='(I2.2)') + $
          STRING(hour, FORMAT='(I2.2)') + $
          STRING(minute, FORMAT='(I2.2)') + $
          STRING(LONG(second), FORMAT='(I2.2)')

  RETURN, stamp
END

; True if token is an integer like 0, 12, +3, -7
FUNCTION is_int_tok_fn, s
  b = BYTE(STRTRIM(s,2))
  IF N_ELEMENTS(b) EQ 0 THEN RETURN, 0
  i = 0L
  ; optional leading sign
  IF (b[0] EQ BYTE('+')) OR (b[0] EQ BYTE('-')) THEN i = 1L
  IF i GE N_ELEMENTS(b) THEN RETURN, 0
  FOR j=i, N_ELEMENTS(b)-1 DO IF (b[j] LT 48B) OR (b[j] GT 57B) THEN RETURN, 0
  RETURN, 1
END

; True if token looks like a float/sci number: 1.23, -3, 4.5e-2, +6E+03, etc.
FUNCTION is_num_tok_fn, s
  t = STRTRIM(s, 2)
  IF t EQ '' THEN RETURN, 0
  b = BYTE(t)
  ; allow digits, ., e/E, +/-, but signs only at start or right after e/E
  seen_digit = 0B
  seen_dot   = 0B
  seen_e     = 0B
  FOR i=0L, N_ELEMENTS(b)-1 DO BEGIN
    c = b[i]
    IF (c GE 48B) AND (c LE 57B) THEN BEGIN
      seen_digit = 1B & CONTINUE
    ENDIF
    IF (c EQ BYTE('.')) THEN BEGIN
      IF seen_dot OR seen_e THEN RETURN, 0
      seen_dot = 1B & CONTINUE
    ENDIF
    IF (c EQ BYTE('e')) OR (c EQ BYTE('E')) THEN BEGIN
      IF seen_e OR ~seen_digit THEN RETURN, 0
      seen_e = 1B
      ; next char may be +/-, but require at least one digit after
      CONTINUE
    ENDIF
    IF (c EQ BYTE('+')) OR (c EQ BYTE('-')) THEN BEGIN
      ; allowed only at start, or immediately after e/E
      IF (i EQ 0) THEN CONTINUE
      IF ((b[i-1] EQ BYTE('e')) OR (b[i-1] EQ BYTE('E'))) THEN CONTINUE
      RETURN, 0
    ENDIF
    ; anything else invalid
    RETURN, 0
  ENDFOR
  RETURN, seen_digit
END

;-------------------------------
; Helper: parse ion tag -> elem, stage, ion_norm
; Returns a structure.
;-------------------------------
FUNCTION clean_ion_fn, tag_in
  COMPILE_OPT idl2
  s = STRLOWCASE(STRTRIM(tag_in, 2))

  p = STRPOS(s, '_')
  IF p LT 0 THEN BEGIN
    elem  = s
    stage = '1'
  ENDIF ELSE BEGIN
    elem  = STRMID(s, 0, p)
    rest  = STRMID(s, p+1)
    b = BYTE(rest)
    idx = WHERE((b GE 48) AND (b LE 57), cnt)
    stage = (cnt GT 0) ? STRJOIN(STRING(b[idx]), '') : '1'
  ENDELSE

  ; keep only a..z in elem
  b = BYTE(elem)
  idx = WHERE((b GE 97) AND (b LE 122), cnt2)
  elem = (cnt2 GT 0) ? STRJOIN(STRING(b[idx]), '') : ''

  ion_norm = elem + '_' + stage
  RETURN, {elem: elem, stage: FIX(LONG(stage)), ion: ion_norm}
END

;-----------------------------------------
; Helper: parse one WGFA data line
; Returns 1 if parsed, 0 otherwise.
;-----------------------------------------
FUNCTION parse_wgfa_line_fn, line, lower, upper, lam, fosc, Aij, descr
  COMPILE_OPT idl2

  s = STRTRIM(line, 2)
  IF s EQ '' THEN RETURN, 0
  IF STRMID(s, 0, 1) EQ ';' THEN RETURN, 0  ; comment lines

  ; Manual whitespace normalize (no STRCOMPRESS/STRTRAN)
  b = BYTE(s)
  out = BYTARR(N_ELEMENTS(b))
  n = 0L
  prev_space = 1B

  FOR i = 0L, N_ELEMENTS(b)-1 DO BEGIN
    c = b[i]
    IF c EQ 9B THEN c = 32B          ; tabs -> space
    IF c LT 32B THEN c = 32B         ; control -> space

    IF c EQ 32B THEN BEGIN
      IF ~prev_space THEN BEGIN
        out[n] = 32B
        n += 1L
        prev_space = 1B
      ENDIF
    ENDIF ELSE BEGIN
      out[n] = c
      n += 1L
      prev_space = 0B
    ENDELSE
  ENDFOR

  IF (n GT 0) AND (out[n-1] EQ 32B) THEN n -= 1L
  IF n LE 0 THEN RETURN, 0
  s2 = STRING(out[0:n-1])

  toks = STRSPLIT(s2, ' ', /EXTRACT)
  IF N_ELEMENTS(toks) LT 6 THEN RETURN, 0

  ; --- validate before converting to avoid IDL errors ---
  IF ~is_int_tok_fn(toks[0]) THEN RETURN, 0
  IF ~is_int_tok_fn(toks[1]) THEN RETURN, 0
  IF ~is_num_tok_fn(toks[2]) THEN RETURN, 0
  IF ~is_num_tok_fn(toks[3]) THEN RETURN, 0
  IF ~is_num_tok_fn(toks[4]) THEN RETURN, 0

  lower = FIX(toks[0])
  upper = FIX(toks[1])
  lam   = DOUBLE(toks[2])
  fosc  = DOUBLE(toks[3])
  Aij   = DOUBLE(toks[4])
  descr = STRJOIN(toks[5:*], ' ')
  RETURN, 1
END


;-----------------------------------------
; Main: find closest wavelength in ion’s WGFA
;-----------------------------------------
FUNCTION chi_find_transition, ion_tag, target_lambda, CHIANTI_ROOT=chi_root, VERBOSE=verbose
  COMPILE_OPT idl2
  IF N_PARAMS() LT 2 THEN MESSAGE, 'Usage: r = chi_find_transition(ion_tag, target_lambda, CHIANTI_ROOT=..., /VERBOSE)'
  IF N_ELEMENTS(verbose) EQ 0 THEN verbose = 0B

  ; --- Resolve CHIANTI root safely ---
  IF N_ELEMENTS(chi_root) EQ 0 THEN chi_root = ''
  IF STRTRIM(chi_root, 2) EQ '' THEN BEGIN
    IF N_ELEMENTS(!xuvtop) NE 0 AND STRTRIM(!xuvtop, 2) NE '' THEN chi_root = !xuvtop
    IF STRTRIM(chi_root, 2) EQ '' THEN BEGIN
      env = GETENV('XUVTOP')
      IF env NE '' THEN chi_root = env
    ENDIF
    IF STRTRIM(chi_root, 2) EQ '' THEN MESSAGE, 'CHIANTI root not set. Pass CHIANTI_ROOT=..., or set !xuvtop or XUVTOP.'
  ENDIF

  ioninfo = clean_ion_fn(ion_tag)
  IF ioninfo.elem EQ '' THEN MESSAGE, 'Could not parse element from ion tag: ' + ion_tag
  ion_norm = ioninfo.ion & elem = ioninfo.elem & stage = ioninfo.stage

  fname = FILEPATH(ion_norm + '.wgfa', ROOT_DIR=chi_root, SUBDIR=[elem, ion_norm])
  IF verbose THEN PRINT, 'Searching WGFA: ', fname

  IF ~FILE_TEST(fname, /REGULAR) THEN RETURN, {found:0B, filepath:fname, ion:ion_norm, element:elem, stage:stage, $
                                               lower:-1L, upper:-1L, wavelength:!VALUES.D_NAN, $
                                               f:!VALUES.D_NAN, A:!VALUES.D_NAN, description:'', delta:!VALUES.D_NAN}

  OPENR, lun, fname, /GET_LUN
  ; best_delta: start with a huge number so first parsed line wins
  best_delta = 1D99   ; avoids needing !VALUES.D_INFINITE
  best_lower = 0L & best_upper = 0L
  best_lam=0D & best_f=0D & best_A=0D & best_descr=''
  line = ''
  WHILE ~EOF(lun) DO BEGIN
    READF, lun, line
    l=0 & u=0 & lam=0D & fosc=0D & Aij=0D & descr=''
    ok = parse_wgfa_line_fn(line, l, u, lam, fosc, Aij, descr)
    IF ok THEN BEGIN
      d = ABS(lam - DOUBLE(target_lambda))
      IF d LT best_delta THEN BEGIN
        best_delta = d & best_lower = l & best_upper = u
        best_lam = lam & best_f = fosc & best_A = Aij & best_descr = descr
      ENDIF
    ENDIF
  ENDWHILE
  FREE_LUN, lun

  found = BYTE(best_delta LT 1D99)
  IF verbose AND found THEN BEGIN
    PRINT, 'Found closest line:'
    PRINT, '  λ = ', STRING(best_lam, FORMAT='(F12.3)'), ' Å   Δλ = ', STRING(best_delta, FORMAT='(F12.3)'), ' Å'
    PRINT, '  lower=', best_lower, '  upper=', best_upper
    PRINT, '  f = ', STRING(best_f, FORMAT='(E12.5)'), '   A = ', STRING(best_A, FORMAT='(E12.5)')
    PRINT, '  ', best_descr
  ENDIF ELSE IF verbose THEN PRINT, 'No parseable data lines found.'

  RETURN, {found:found, filepath:fname, ion:ion_norm, element:elem, stage:stage, $
           lower:best_lower, upper:best_upper, wavelength:best_lam, f:best_f, A:best_A, $
           description:best_descr, delta:best_delta}
END
