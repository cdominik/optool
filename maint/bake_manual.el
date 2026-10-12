  (eval-after-load "ox-latex"
    '(add-to-list 'org-latex-classes
                  '("koma-article" "\\documentclass{scrartcl}"
                    ("\\section{%s}" . "\\section*{%s}")
                    ("\\subsection{%s}" . "\\subsection*{%s}")
                    ("\\subsubsection{%s}" . "\\subsubsection*{%s}")
                    ("\\paragraph{%s}" . "\\paragraph*{%s}")
                    ("\\subparagraph{%s}" . "\\subparagraph*{%s}"))))
  ;; The KOMA-Script classes define \captionof themselves, loading
  ;; capt-of on top of that is an error.
  (eval-after-load "ox-latex"
    '(setq org-latex-default-packages-alist
           (seq-remove (lambda (p) (equal (cadr p) "capt-of"))
                       org-latex-default-packages-alist)))
