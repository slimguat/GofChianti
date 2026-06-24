; I can't specify th temperature array it is always the same so I will remove all the T input
pro compute_gofnt
  compile_opt idl2
  on_error, 2 ; stop on error

  required = [ $
    'datetime_stamp', $
    'is_int_tok_fn', $
    'is_num_tok_fn', $
    'clean_ion_fn', $
    'parse_wgfa_line_fn', $
    'chi_find_transition' $
    ]

  for k = 0l, n_elements(required) - 1 do begin
    catch, err
    if err ne 0 then begin
      catch, /cancel
      message, 'Missing or failed routine: ' + required[k] + $
        ' (ensure the .pro is on !PATH and compiles).'
    endif

    ; Try to compile from source quietly; if it fails, CATCH handles it
    resolve_routine, required[k], /compile, /quiet
    catch, /cancel
  endfor

  ; ---------------- User parameters ----------------
  VERBOSE = 1b
  Hold_for_debug = 0b
  if VERBOSE then print, 'Starting GOFNT computations...'
  if VERBOSE then print, 'VERBOSE mode is', VERBOSE, 'Hold_for_debug mode is', Hold_for_debug
  OUTDIR = './gofnt_remove_later/'
  if ~file_test(OUTDIR, /dir) then file_mkdir, OUTDIR
  if VERBOSE then print, 'GOFNT output directory: ', OUTDIR

  redo_calculation = 0b ; set 1B to recompute even if file exists

  ; lines = ['N_4 765.152']
  ; lines = [ $
  ; 'N_4 765.152', $
  ; 'S_4 750.221', $
  ; 'S_5 786.468', $
  ; 'N_3 991.511', $
  ; 'N_3 991.577'  $
  ; ]
  ; SPICE spectral lines
  ; -----------------------------------------
  ; Composition and dynamics lines combined
  ; -----------------------------------------

  lines = [ $
    'H_1 1025.723' $
    ; ; EIS lines
    ; 'Fe_12 195.119', $
    ; 'Si_10 258.374', $
    ; 'S_10 264.230', $
    ; 'Fe_13 202.044', $
    ; 'Fe_14 270.52', $
    ; ; 'Ar_8 700.240', $
    ; ; Composition lines
    ; 'S_4 750.221', $
    ; 'N_4 765.152', $
    ; 'S_5 786.468', $
    ; 'N_3 991.511', $
    ; 'N_3 991.577', $
    ; 'N_3 989.799', $

    ; ; Density lines
    ; 'O_5 758.677', $
    ; 'O_5 759.442', $
    ; 'O_5 760.227', $
    ; 'O_5 760.446', $
    ; 'O_5 761.128', $
    ; 'O_5 762.004', $
    ; 'O_5 764.561', $

    ; ; Temperature lines
    ; 'O_3 702.838', $
    ; 'O_3 702.896', $
    ; 'O_3 702.337', $
    ; 'O_3 702.900', $
    ; 'O_3 703.851', $
    ; 'O_3 703.855', $
    ; 'O_4 787.710', $

    ; ; Dynamics lines
    ; 'Mg_9 706.060', $
    ; 'Mg_9 749.552', $
    ; 'S_4 748.393', $
    ; 'Mg_8 769.355', $
    ; 'Ne_8 770.428', $
    ; 'Mg_8 772.260', $
    ; 'Ne_8 780.385', $
    ; 'Mg_8 782.362', $
    ; 'Fe_18 974.858', $
    ; 'C_3 977.020', $

    ; ; Extra weak lines
    ; 'O_5 760.674', $
    ; 'O_5 761.821', $
    ; 'O_5 761.815', $
    ; 'O_5 763.662', $
    ; 'O_5 763.657', $
    ; 'O_5 763.644', $
    ; 'O_5 763.942', $
    ; 'O_5 763.592', $
    ; 'O_5 764.584', $
    ; 'O_5 764.935', $
    ; 'O_4 790.114', $
    ; 'O_4 790.201', $
    ; 'O_4 790.201', $
    ; ; ; Weird lines
    ; 'H_1 972.537', $ ; (two transitions: 1→11 and 1→12)
    ; ; 'He_2 303.786', $
    ; ; 'He_2 303.785' $
    ; ; Extra lines
    ; 'ar_8 700.240', $
    ; 'fe_3 985.852', $
    ; 'ne_6 992.683', $
    ; 'Fe_3 994.258', $
    ; 'na_6 988.709' $
    ]

  nN = 10l ; (7-1)*33 + 1 -> heads at 1e7..1e13
  Nmin = 1.0e7
  Nmax = 1.0e13
  ; nN = 66l ; (7-1)*33 + 1 -> heads at 1e7..1es13
  ; Nmin = 1.0e1
  ; Nmax = 1.0e14

  ; nN = 2L   ; (7-1)*33 + 1 -> heads at 1e7..1e13
  ; Nmin = 1.0e9  & Nmax = 1.0e10

  if VERBOSE then begin
    print, 'Density    grid: log10(n) ', Nmin, ' -> ', Nmax, ' (cm^-3), nN=', nN
  endif

  ; Ensure CHIANTI env before using !xuvtop
  if (n_elements(!xuvtop) eq 0) or (strlen(!xuvtop) eq 0) then message, 'ERROR: !XUVTOP not set.'

  IONEQ_NAME = './ioneq_tmp/ioneq' + datetime_stamp() + '/'
  if ~file_test(IONEQ_NAME, /dir) then file_mkdir, IONEQ_NAME
  if VERBOSE then print, 'Using IONEQ_NAME directory: ', IONEQ_NAME

  ABUND_NAME = !xuvtop + '/abundance/sun_photospheric_2021_asplund.abund'
  ; Store the BARE contribution function (no abundance). The Python GofChianti
  ; package multiplies by the chosen abundance at read time. Set to 0b only if
  ; you specifically want abundance baked into the .dat output.
  NOABUND = 1b
  ; ---------------- End user parameters ------------

  ; ---- CHIANTI version ----
  filename = !xuvtop + '/VERSION'
  openr, lun_v, filename, /get_lun
  first_line = ''
  readf, lun_v, first_line
  free_lun, lun_v
  chianti_version = strtrim(first_line, 2)
  if VERBOSE then print, 'CHIANTI version: ', chianti_version

  ; ---- Grids ----
  logN = findgen(nN) / (nN - 1) * (alog10(Nmax) - alog10(Nmin)) + alog10(Nmin)
  dens = 10.0 ^ logN
  nN = n_elements(dens)

  if Hold_for_debug then begin
    ; prompt input once clicked enter it will start processing
    print, 'Press 0+Enter to start GOFNT calculations...'
    read, dummy_input
  endif

  ; ---- Loop over requested lines ----
  for i = 0l, n_elements(lines) - 1 do begin
    parts = strsplit(lines[i], ' ', /extract)
    ion_tag = parts[0]
    wvl = float(parts[1])

    ; Find transition
    transition = chi_find_transition(ion_tag, wvl, verbose = VERBOSE)

    ; if not found skip to next
    if transition.found eq 0b then begin
      print, '  Transition not found for ', ion_tag, ' at ', string(wvl, format = '(F7.3)'), ' Å. Skipping.'
      continue
    endif

    ion_tag = transition.ion

    print, 'Processing line: ', ion_tag, ' at ', string(wvl, format = '(F7.3)'), ' Å  => λ_found = ', $
      string(transition.wavelength, format = '(F7.3)'), ' Å'

    ; Output filename
    wvl_str = strcompress(string(transition.wavelength, format = '(F8.3)'), /remove_all)
    save_file = OUTDIR + ion_tag + '_' + wvl_str + '_gofnt_v-' + chianti_version + '.dat'

    ; Skip if exists and not redoing
    if file_test(save_file) and (redo_calculation eq 0b) then begin
      print, '  Exists, skipping: ', save_file
      continue
    endif

    if VERBOSE then begin
      print, '  Ion tag: ', ion_tag
      print, '  Wmin = ', string(transition.wavelength - 5, format = '(F7.3)'), ' Å'
      print, '  Wmax = ', string(transition.wavelength + 5, format = '(F7.3)'), ' Å'
      print, '  NOABUND = ', NOABUND
      print, '  ABUND_NAME = ', ABUND_NAME
      print, '  IONEQ_NAME_FOLDER = ', IONEQ_NAME
      print, '  upper_levels = ', transition.upper
      print, '  lower_levels = ', transition.lower
    endif

    ; Compute g(T) for each density
    if Hold_for_debug then begin
      print, 'Press 0+Enter to start GOFNT calculations for this line...'
      read, dummy_input
    endif
    for j = 0l, nN - 1 do begin
      ; print j value
      if VERBOSE then print, '  Computing GOFNT for density ', string(dens[j], format = '(E12.5)'), ' cm^-3  (', j + 1, ' of ', nN, ')'

      timestr = datetime_stamp()
      normal_date = strmid(timestr, 0, 4) + '-' + strmid(timestr, 4, 2) + '-' + $
        strmid(timestr, 6, 2) + ' ' + strmid(timestr, 9, 2) + ':' + $
        strmid(timestr, 11, 2) + ':' + strmid(timestr, 13, 2)
      sub_IONEQ_NAME = IONEQ_NAME + timestr + '.ioneq'
      ; Slight window around the found wavelength
      gofnt, ion_tag, transition.wavelength - 5, transition.wavelength + 5, temp, g, desc, $
        density = dens[j], noabund = NOABUND, $
        abund_name = ABUND_NAME, ioneq_name = sub_IONEQ_NAME, $
        upper_levels = transition.upper, lower_levels = transition.lower

      if j eq 0 then begin
        ; open a new file and dump the header
        openw, lun, save_file, /get_lun
        if VERBOSE then print, '  Saving to file: ', save_file

        ; Header
        ; NOTE: each line is built as a SINGLE string and printed with format
        ; '(A)' so that (a) long values (e.g. the abundance path) are never
        ; wrapped at IDL's default 80-column limit, and (b) the wavelength is
        wl_str = strtrim(string(transition.wavelength, format = '(F12.3)'), 2)
        printf, lun, format = '(A)', '# GOFNT data for line: ' + ion_tag + ' at ' + wl_str + ' Å'
        printf, lun, format = '(A)', '# Generated on: ' + normal_date
        printf, lun, format = '(A)', '# CHIANTI version: ' + chianti_version
        printf, lun, format = '(A)', '# nN: ' + strtrim(string(nN), 2)
        printf, lun, format = '(A)', '# f: ' + strtrim(string(transition.f, format = '(E12.5)'), 2)
        printf, lun, format = '(A)', '# A: ' + strtrim(string(transition.a, format = '(E12.5)'), 2)
        printf, lun, format = '(A)', '# description: ' + transition.description
        printf, lun, format = '(A)', '# LOWER LEVEL: ' + strtrim(string(transition.lower), 2)
        printf, lun, format = '(A)', '# UPPER LEVEL: ' + strtrim(string(transition.upper), 2)
        printf, lun, format = '(A)', '# minN: ' + strtrim(string(Nmin, format = '(E12.5)'), 2) + ' cm^-3'
        printf, lun, format = '(A)', '# maxN: ' + strtrim(string(Nmax, format = '(E12.5)'), 2) + ' cm^-3'
        ; Only the basename (no path): users do not care about the SSW path and
        ; the full path is what used to trigger the 80-column line wrap.
        printf, lun, format = '(A)', '# Abundance file: ' + ABUND_NAME
        printf, lun, format = '(A)', '# Abundance multiplication: ' + strtrim(string(fix(~NOABUND)), 2)
        printf, lun, format = '(A)', '# Rows and units: Density (cm^-3), Temperature (K), G(T,n) (erg cm^3 s^-1 sr^-1)'
      endif

      ; Build single-line strings before printing
      d_str = string(dens[j], format = '(E12.5)')
      t_str = strjoin(string(temp, format = '(E12.5)'), ' ')
      g_str = strjoin(string(g, format = '(E12.5)'), ' ')
      ; Now print them line by line
      printf, lun, 'De: ', d_str
      printf, lun, t_str
      printf, lun, g_str

      ; gofnt, ion_tag, transition.wavelength - 5, transition.wavelength + 5, T2, g, desc, $
      ; DENSITY = dens[j], $
      ; ABUND_NAME = ABUND_NAME,LOGT0= log_T, $
      ; upper_levels = transition.upper, lower_levels = transition.lower;,VERBOSE= VERBOSE
    endfor

    free_lun, lun
    if VERBOSE then print, '  GOFNT calculation completed and saved for ', ion_tag, ' at ', $
      string(wvl, format = '(F7.3)'), ' Å'
  endfor
end
