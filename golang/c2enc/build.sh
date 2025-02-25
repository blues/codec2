#!/bin/bash

go fmt
go build
go vet
staticcheck
