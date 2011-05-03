##-- CUSTOM Cufflinks RULES --##

##-- install-exec-hook
##   execute this custom rule after the global install-exec rule
##
install-exec-hook: install-bin-scripts-hook

##-- uninstall-hook
##   execute this custom rule after the global uninstall rule
##
uninstall-hook: uninstall-bin-scripts-hook

##-- install-bin-scripts-hook
##   takes list of scripts from $(dist_bin_SCRIPTS) $(bin_SCRIPTS) and:
##   1) locates the installed script
##   2) builds them with the below sed commands to assure the shebang and
##   include directives point to the correct locations, and that variable
##   substitution is performed as expected
##   3) strips their suffix (e.g. "script.sh" becomes "script")
##
install-bin-scripts-hook: $(dist_bin_SCRIPTS) $(bin_SCRIPTS)
@$(NORMAL_INSTALL)
@list='$(dist_bin_SCRIPTS) $(bin_SCRIPTS)'; for p in $$list; do \
f=`echo "$$p" | sed 's|^.*/||;$(transform)'`; \
b=`echo "$$p" | sed 's|^.*/||;s|\.[^.]*$$||;$(transform)'`; \
if test -f "$(DESTDIR)$(bindir)/$$f"  &&  test "$$f" != "$$b"; then \
case "$$p" in \
*.py) \
echo " configuring python '$$b'"; \
echo '#!$(PYTHON)' > "$(DESTDIR)$(bindir)/$$b"; \
sed -e '1 {s|^#!.*$$||;}' \
"$(DESTDIR)$(bindir)/$$f" >> "$(DESTDIR)$(bindir)/$$b" \
|| exit 1; \
;; \
*.sh) \
echo " configuring shell '$$b', using #!$(TOPHAT_SHELL)"; \
echo '#!$(TOPHAT_SHELL)' > "$(DESTDIR)$(bindir)/$$b"; \
sed -e '1 {s|^#!.*$$||;}' \
-e 's|BINDIR[[:space:]]*=.*|BINDIR=$(bindir)|' \
"$(DESTDIR)$(bindir)/$$f" >> "$(DESTDIR)$(bindir)/$$b" \
|| exit 1; \
;; \
*.pl) \
echo " configuring perl '$$b'"; \
echo '#!$(PERL)' > "$(DESTDIR)$(bindir)/$$b"; \
sed -e '1 {s|^#!.*$$||;}' \
"$(DESTDIR)$(bindir)/$$f" >> "$(DESTDIR)$(bindir)/$$b" \
|| exit 1; \
;; \
*) \
echo " configuring '$$b'"; \
cp "$(DESTDIR)$(bindir)/$$f" "$(DESTDIR)$(bindir)/$$b" \
;; \
esac; \
mv -f "$(DESTDIR)$(bindir)/$$b" "$(DESTDIR)$(bindir)/$$f" || exit 1; \
$(INSTALL_SCRIPT) "$(DESTDIR)$(bindir)/$$f" "$(DESTDIR)$(bindir)/$$b" \
|| exit 1; \
rm -f "$(DESTDIR)$(bindir)/$$f"; \
else :; fi; \
done


##-- uninstall-bin-scripts-hook
##   takes list of scripts from $(dist_bin_SCRIPTS) $(bin_SCRIPTS) and removes
##   their suffix-stripped versions.
##
uninstall-bin-scripts-hook:
@$(NORMAL_UNINSTALL)
@list='$(dist_bin_SCRIPTS) $(bin_SCRIPTS)'; for p in $$list; do \
b=`echo "$$p" | sed 's|^.*/||;s|.[^.]*$$||;$(transform)'`; \
echo " rm -f '$(DESTDIR)$(bindir)/$$b'"; \
rm -f "$(DESTDIR)$(bindir)/$$b"; \
done


##-- END OF MAKEFILE --##