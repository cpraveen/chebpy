html:
	export BASE_URL=/chebpy
	myst build --html --execute

upload:
	ghp-import -n -p -f _build/html

clean:
	rm -rf _build
