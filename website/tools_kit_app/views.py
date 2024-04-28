from django.shortcuts import render

# Create your views here.
#request handler

from django.https import HttpResponse

def say_hello(request):
    return HttpResponse('Hello World')