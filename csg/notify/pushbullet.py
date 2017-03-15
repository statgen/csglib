import os, sys, psutil, time, requests
from subprocess import Popen

def send_message(email,title,message,token):
  """
  Send a notification message using Pushbullet. Requires API token.

  Example: 
    send_message("Job finished","My job finally finished","FujY432FnBAhrueAFD312141")

  :param title: Title of notification message
  :param message: Body of the notification message
  :param token: This is your Pushbullet API token. You must register for a Pushbullet account to obtain one. 
  """

  headers = {
    "Access-Token" : token
  }

  payload = {
    "type" : "note",
    "body" : message,
    "title" : title,
    "email" : email
  }

  req = requests.post(
    "https://api.pushbullet.com/v2/pushes",
    headers = headers,
    json = payload
  )

  return req

