html:
	export BASE_URL=/chebpy
	jupyter book build --html

upload:
	ghp-import -n -p -f _build/html

clean:
	rm -rf _build
