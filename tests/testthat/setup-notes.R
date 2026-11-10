## coloc() and wakefield_bf() note a plain list once per session. The tests run
## as a session that has had the note, so the many that pass plain lists stay
## quiet whatever order the files run in; test-notes.R forgets it where it
## tests it. A setup file, not a helper, so that devtools::load_all() leaves
## the note to be given
.coloc_state$noted <- union(.coloc_state$noted, "plain_list")
